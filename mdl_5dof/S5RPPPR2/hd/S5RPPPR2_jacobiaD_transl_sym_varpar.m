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
% Datum: 2019-10-24 10:39
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
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
	% StartTime: 2019-10-24 10:39:20
	% EndTime: 2019-10-24 10:39:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:39:20
	% EndTime: 2019-10-24 10:39:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; (r_i_i_C(1) * t5 + r_i_i_C(2) * t6) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t6 + r_i_i_C(2) * t5) * qJD(1), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:39:20
	% EndTime: 2019-10-24 10:39:20
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (11->8), mult. (28->12), div. (0->0), fcn. (18->4), ass. (0->5)
	t37 = -r_i_i_C(3) - qJ(2);
	t36 = r_i_i_C(1) * cos(pkin(7)) - r_i_i_C(2) * sin(pkin(7)) + pkin(1);
	t35 = cos(qJ(1));
	t34 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; -t34 * qJD(2) + (t36 * t34 + t37 * t35) * qJD(1), -qJD(1) * t34, 0, 0, 0; t35 * qJD(2) + (t37 * t34 - t36 * t35) * qJD(1), qJD(1) * t35, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:39:20
	% EndTime: 2019-10-24 10:39:20
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (23->16), mult. (62->21), div. (0->0), fcn. (48->6), ass. (0->11)
	t75 = sin(pkin(8));
	t76 = sin(pkin(7));
	t77 = cos(pkin(8));
	t87 = (r_i_i_C(3) + qJ(3)) * t76 + (r_i_i_C(1) * t77 - r_i_i_C(2) * t75 + pkin(2)) * cos(pkin(7)) + pkin(1);
	t79 = sin(qJ(1));
	t85 = qJD(1) * t79;
	t80 = cos(qJ(1));
	t84 = qJD(1) * t80;
	t83 = t76 * qJD(3);
	t81 = -t75 * r_i_i_C(1) - t77 * r_i_i_C(2) - qJ(2);
	t1 = [0, 0, 0, 0, 0; -t80 * t83 - t79 * qJD(2) + (t87 * t79 + t81 * t80) * qJD(1), -t85, -t76 * t84, 0, 0; -t79 * t83 + t80 * qJD(2) + (t81 * t79 - t87 * t80) * qJD(1), t84, -t76 * t85, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:39:20
	% EndTime: 2019-10-24 10:39:20
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (43->29), mult. (126->39), div. (0->0), fcn. (110->8), ass. (0->23)
	t110 = sin(pkin(9));
	t112 = sin(pkin(7));
	t113 = cos(pkin(9));
	t115 = cos(pkin(7));
	t131 = pkin(2) * t115 + pkin(1) + (r_i_i_C(1) * t110 + r_i_i_C(2) * t113 + qJ(3)) * t112;
	t130 = r_i_i_C(3) + qJ(4);
	t111 = sin(pkin(8));
	t116 = sin(qJ(1));
	t129 = t111 * t116;
	t117 = cos(qJ(1));
	t128 = t111 * t117;
	t114 = cos(pkin(8));
	t127 = t116 * t114;
	t126 = t117 * t114;
	t125 = qJD(1) * t116;
	t124 = qJD(1) * t117;
	t123 = t112 * qJD(3);
	t122 = t115 * t128;
	t119 = t115 * t129 + t126;
	t118 = qJD(1) * (t113 * r_i_i_C(1) - t110 * r_i_i_C(2) + pkin(3));
	t107 = qJD(1) * t122 - t114 * t125;
	t105 = t119 * qJD(1);
	t1 = [0, 0, 0, 0, 0; -(t122 - t127) * qJD(4) - t117 * t123 - t116 * qJD(2) + t130 * t105 + (t115 * t127 - t128) * t118 + (-t117 * qJ(2) + t131 * t116) * qJD(1), -t125, -t112 * t124, -t107, 0; -t119 * qJD(4) - t116 * t123 + t117 * qJD(2) - t130 * t107 + (-t115 * t126 - t129) * t118 + (-t116 * qJ(2) - t131 * t117) * qJD(1), t124, -t112 * t125, -t105, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:39:20
	% EndTime: 2019-10-24 10:39:20
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (128->59), mult. (392->96), div. (0->0), fcn. (392->10), ass. (0->44)
	t232 = cos(pkin(7));
	t228 = sin(pkin(8));
	t236 = cos(qJ(1));
	t246 = t236 * t228;
	t231 = cos(pkin(8));
	t234 = sin(qJ(1));
	t247 = t234 * t231;
	t237 = t232 * t247 - t246;
	t218 = t237 * qJD(1);
	t227 = sin(pkin(9));
	t230 = cos(pkin(9));
	t229 = sin(pkin(7));
	t244 = qJD(1) * t234;
	t240 = t229 * t244;
	t210 = t218 * t230 + t227 * t240;
	t245 = t236 * t231;
	t248 = t234 * t228;
	t222 = t232 * t248 + t245;
	t217 = t222 * qJD(1);
	t233 = sin(qJ(5));
	t235 = cos(qJ(5));
	t259 = t210 * t233 - t217 * t235;
	t258 = t210 * t235 + t217 * t233;
	t224 = t232 * t245 + t248;
	t249 = t229 * t227;
	t215 = t224 * t230 + t236 * t249;
	t241 = t232 * t246;
	t223 = t241 - t247;
	t257 = t215 * t233 - t223 * t235;
	t256 = t215 * t235 + t223 * t233;
	t255 = r_i_i_C(3) + pkin(6);
	t250 = t228 * t229;
	t243 = qJD(1) * t236;
	t242 = t229 * qJD(3);
	t239 = t229 * t243;
	t238 = pkin(2) * t232 + qJ(3) * t229 + pkin(1);
	t221 = t229 * t231 * t230 - t232 * t227;
	t220 = t224 * qJD(1);
	t219 = qJD(1) * t241 - t231 * t244;
	t214 = -t237 * t230 - t234 * t249;
	t213 = -t220 * t230 - t227 * t239;
	t208 = t213 * t235 - t219 * t233 + (-t214 * t233 - t222 * t235) * qJD(5);
	t207 = -t213 * t233 - t219 * t235 + (-t214 * t235 + t222 * t233) * qJD(5);
	t1 = [0, 0, 0, 0, ((-t221 * t235 - t233 * t250) * r_i_i_C(1) + (t221 * t233 - t235 * t250) * r_i_i_C(2)) * qJD(5); t258 * r_i_i_C(1) - t259 * r_i_i_C(2) + t210 * pkin(4) + t218 * pkin(3) + t217 * qJ(4) - t223 * qJD(4) - t236 * t242 - t234 * qJD(2) + t255 * (t218 * t227 - t230 * t240) + (t257 * r_i_i_C(1) + t256 * r_i_i_C(2)) * qJD(5) + (-t236 * qJ(2) + t238 * t234) * qJD(1), -t244, -t239, -t219, t207 * r_i_i_C(1) - t208 * r_i_i_C(2); -t234 * t242 - t220 * pkin(3) + t213 * pkin(4) + t208 * r_i_i_C(1) + t207 * r_i_i_C(2) - t219 * qJ(4) + t236 * qJD(2) - t222 * qJD(4) + t255 * (-t220 * t227 + t230 * t239) + (-qJ(2) * t234 - t238 * t236) * qJD(1), t243, -t240, -t217, t259 * r_i_i_C(1) + t258 * r_i_i_C(2) + (-t256 * r_i_i_C(1) + t257 * r_i_i_C(2)) * qJD(5);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end