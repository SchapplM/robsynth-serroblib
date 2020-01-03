% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPPPR2
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPPPR2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPPPR2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:23:58
	% EndTime: 2020-01-03 11:23:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:23:58
	% EndTime: 2020-01-03 11:23:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; (-r_i_i_C(1) * t5 - r_i_i_C(2) * t6) * qJD(1), 0, 0, 0, 0; (r_i_i_C(1) * t6 - r_i_i_C(2) * t5) * qJD(1), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:23:58
	% EndTime: 2020-01-03 11:23:58
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (11->8), mult. (28->12), div. (0->0), fcn. (18->4), ass. (0->5)
	t37 = r_i_i_C(3) + qJ(2);
	t36 = r_i_i_C(1) * cos(pkin(7)) - r_i_i_C(2) * sin(pkin(7)) + pkin(1);
	t35 = cos(qJ(1));
	t34 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; t34 * qJD(2) + (-t36 * t34 + t37 * t35) * qJD(1), qJD(1) * t34, 0, 0, 0; -t35 * qJD(2) + (t37 * t34 + t36 * t35) * qJD(1), -qJD(1) * t35, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:23:58
	% EndTime: 2020-01-03 11:23:58
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (21->14), mult. (62->21), div. (0->0), fcn. (48->6), ass. (0->11)
	t71 = sin(pkin(8));
	t72 = sin(pkin(7));
	t73 = cos(pkin(8));
	t83 = (r_i_i_C(3) + qJ(3)) * t72 + (r_i_i_C(1) * t73 - r_i_i_C(2) * t71 + pkin(2)) * cos(pkin(7)) + pkin(1);
	t75 = sin(qJ(1));
	t81 = qJD(1) * t75;
	t76 = cos(qJ(1));
	t80 = qJD(1) * t76;
	t79 = t72 * qJD(3);
	t77 = t71 * r_i_i_C(1) + t73 * r_i_i_C(2) + qJ(2);
	t1 = [0, 0, 0, 0, 0; t76 * t79 + t75 * qJD(2) + (-t83 * t75 + t77 * t76) * qJD(1), t81, t72 * t80, 0, 0; t75 * t79 - t76 * qJD(2) + (t77 * t75 + t83 * t76) * qJD(1), -t80, t72 * t81, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:23:58
	% EndTime: 2020-01-03 11:23:58
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (41->24), mult. (126->39), div. (0->0), fcn. (110->8), ass. (0->22)
	t116 = sin(pkin(9));
	t118 = sin(pkin(7));
	t119 = cos(pkin(9));
	t121 = cos(pkin(7));
	t136 = pkin(2) * t121 + pkin(1) + (r_i_i_C(1) * t116 + r_i_i_C(2) * t119 + qJ(3)) * t118;
	t135 = r_i_i_C(3) + qJ(4);
	t117 = sin(pkin(8));
	t122 = sin(qJ(1));
	t134 = t117 * t122;
	t123 = cos(qJ(1));
	t133 = t121 * t123;
	t120 = cos(pkin(8));
	t132 = t122 * t120;
	t131 = qJD(1) * t122;
	t130 = qJD(1) * t123;
	t129 = qJD(3) * t118;
	t126 = t117 * t133 - t132;
	t125 = t120 * t123 + t121 * t134;
	t124 = qJD(1) * (r_i_i_C(1) * t119 - r_i_i_C(2) * t116 + pkin(3));
	t114 = t126 * qJD(1);
	t112 = t125 * qJD(1);
	t1 = [0, 0, 0, 0, 0; t126 * qJD(4) + t123 * t129 + t122 * qJD(2) - t135 * t112 + (t117 * t123 - t121 * t132) * t124 + (t123 * qJ(2) - t136 * t122) * qJD(1), t131, t118 * t130, t114, 0; t125 * qJD(4) + t122 * t129 - t123 * qJD(2) + t135 * t114 + (t120 * t133 + t134) * t124 + (t122 * qJ(2) + t136 * t123) * qJD(1), -t130, t118 * t131, t112, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:23:58
	% EndTime: 2020-01-03 11:23:59
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (126->55), mult. (392->96), div. (0->0), fcn. (392->10), ass. (0->43)
	t231 = sin(pkin(8));
	t235 = cos(pkin(7));
	t234 = cos(pkin(8));
	t237 = sin(qJ(1));
	t249 = t237 * t234;
	t243 = t235 * t249;
	t239 = cos(qJ(1));
	t245 = qJD(1) * t239;
	t221 = qJD(1) * t243 - t231 * t245;
	t230 = sin(pkin(9));
	t233 = cos(pkin(9));
	t232 = sin(pkin(7));
	t246 = qJD(1) * t237;
	t242 = t232 * t246;
	t212 = t221 * t233 + t230 * t242;
	t247 = t239 * t234;
	t250 = t237 * t231;
	t225 = t235 * t250 + t247;
	t220 = t225 * qJD(1);
	t236 = sin(qJ(5));
	t238 = cos(qJ(5));
	t261 = t212 * t236 - t220 * t238;
	t260 = -t212 * t238 - t220 * t236;
	t227 = t235 * t247 + t250;
	t251 = t232 * t230;
	t218 = t227 * t233 + t239 * t251;
	t248 = t239 * t231;
	t226 = t235 * t248 - t249;
	t259 = -t218 * t236 + t226 * t238;
	t258 = t218 * t238 + t226 * t236;
	t257 = r_i_i_C(3) + pkin(6);
	t252 = t231 * t232;
	t244 = t232 * qJD(3);
	t241 = t232 * t245;
	t240 = pkin(2) * t235 + qJ(3) * t232 + pkin(1);
	t224 = t232 * t234 * t233 - t235 * t230;
	t223 = t227 * qJD(1);
	t222 = t226 * qJD(1);
	t217 = (t243 - t248) * t233 + t237 * t251;
	t216 = t223 * t233 + t230 * t241;
	t211 = t216 * t238 + t222 * t236 + (-t217 * t236 + t225 * t238) * qJD(5);
	t210 = -t216 * t236 + t222 * t238 + (-t217 * t238 - t225 * t236) * qJD(5);
	t1 = [0, 0, 0, 0, ((-t224 * t238 - t236 * t252) * r_i_i_C(1) + (t224 * t236 - t238 * t252) * r_i_i_C(2)) * qJD(5); t260 * r_i_i_C(1) + t261 * r_i_i_C(2) - t212 * pkin(4) - t221 * pkin(3) - t220 * qJ(4) + t226 * qJD(4) + t239 * t244 + t237 * qJD(2) + t257 * (-t221 * t230 + t233 * t242) + (t259 * r_i_i_C(1) - t258 * r_i_i_C(2)) * qJD(5) + (t239 * qJ(2) - t240 * t237) * qJD(1), t246, t241, t222, t210 * r_i_i_C(1) - t211 * r_i_i_C(2); t237 * t244 + t223 * pkin(3) + t216 * pkin(4) + t211 * r_i_i_C(1) + t210 * r_i_i_C(2) + t222 * qJ(4) - t239 * qJD(2) + t225 * qJD(4) + t257 * (t223 * t230 - t233 * t241) + (qJ(2) * t237 + t240 * t239) * qJD(1), -t245, t242, t220, -t261 * r_i_i_C(1) + t260 * r_i_i_C(2) + (t258 * r_i_i_C(1) + t259 * r_i_i_C(2)) * qJD(5);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end