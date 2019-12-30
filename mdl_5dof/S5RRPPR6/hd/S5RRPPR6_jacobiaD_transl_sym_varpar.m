% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPPR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRPPR6_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR6_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPPR6_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR6_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:19:54
	% EndTime: 2019-12-29 18:19:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:19:49
	% EndTime: 2019-12-29 18:19:49
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:19:49
	% EndTime: 2019-12-29 18:19:49
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (19->15), mult. (64->29), div. (0->0), fcn. (40->4), ass. (0->13)
	t28 = pkin(6) + r_i_i_C(3);
	t18 = sin(qJ(1));
	t27 = qJD(1) * t18;
	t20 = cos(qJ(1));
	t26 = qJD(1) * t20;
	t25 = qJD(2) * t18;
	t24 = qJD(2) * t20;
	t17 = sin(qJ(2));
	t19 = cos(qJ(2));
	t23 = r_i_i_C(1) * t17 + r_i_i_C(2) * t19;
	t22 = -r_i_i_C(1) * t19 + r_i_i_C(2) * t17 - pkin(1);
	t21 = t23 * qJD(2);
	t1 = [t23 * t25 + (-t18 * t28 + t22 * t20) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0, 0; -t20 * t21 + (t22 * t18 + t20 * t28) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0, 0; 0, -t21, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:19:49
	% EndTime: 2019-12-29 18:19:49
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (50->16), mult. (90->25), div. (0->0), fcn. (59->6), ass. (0->14)
	t27 = qJ(2) + pkin(8);
	t25 = sin(t27);
	t26 = cos(t27);
	t42 = -r_i_i_C(1) * t26 + r_i_i_C(2) * t25 - cos(qJ(2)) * pkin(2);
	t40 = r_i_i_C(3) + qJ(3) + pkin(6);
	t32 = cos(qJ(1));
	t39 = qJD(1) * t32;
	t37 = -pkin(1) + t42;
	t36 = sin(qJ(2)) * pkin(2) + r_i_i_C(1) * t25 + r_i_i_C(2) * t26;
	t30 = sin(qJ(1));
	t35 = t36 * t30;
	t34 = qJD(2) * t42;
	t33 = t36 * qJD(2);
	t1 = [t32 * qJD(3) + qJD(2) * t35 + (-t40 * t30 + t37 * t32) * qJD(1), qJD(1) * t35 + t32 * t34, t39, 0, 0; t30 * qJD(3) - t32 * t33 + (t37 * t30 + t40 * t32) * qJD(1), t30 * t34 - t36 * t39, qJD(1) * t30, 0, 0; 0, -t33, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:19:55
	% EndTime: 2019-12-29 18:19:56
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (123->22), mult. (204->32), div. (0->0), fcn. (152->8), ass. (0->19)
	t176 = sin(pkin(9));
	t177 = cos(pkin(9));
	t175 = qJ(2) + pkin(8);
	t173 = sin(t175);
	t174 = cos(t175);
	t190 = r_i_i_C(1) * t177 - r_i_i_C(2) * t176 + pkin(3);
	t195 = r_i_i_C(3) + qJ(4);
	t186 = t190 * t173 - t195 * t174 + sin(qJ(2)) * pkin(2);
	t183 = -t186 * qJD(2) + t173 * qJD(4);
	t200 = t183 + (t176 * r_i_i_C(1) + t177 * r_i_i_C(2) + pkin(6) + qJ(3)) * qJD(1);
	t199 = -t195 * t173 - t190 * t174 - cos(qJ(2)) * pkin(2);
	t180 = sin(qJ(1));
	t194 = qJD(1) * t180;
	t182 = cos(qJ(1));
	t193 = qJD(1) * t182;
	t192 = qJD(2) * t174;
	t185 = qJD(3) + (-pkin(1) + t199) * qJD(1);
	t184 = t199 * qJD(2) + qJD(4) * t174;
	t1 = [-t200 * t180 + t185 * t182, t184 * t182 + t186 * t194, t193, -t173 * t194 + t182 * t192, 0; t185 * t180 + t200 * t182, t184 * t180 - t186 * t193, t194, t173 * t193 + t180 * t192, 0; 0, t183, 0, qJD(2) * t173, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:19:51
	% EndTime: 2019-12-29 18:19:51
	% DurationCPUTime: 0.41s
	% Computational Cost: add. (264->46), mult. (329->72), div. (0->0), fcn. (261->10), ass. (0->37)
	t215 = cos(pkin(9)) * pkin(4) + pkin(3);
	t222 = qJ(2) + pkin(8);
	t218 = sin(t222);
	t220 = cos(t222);
	t253 = r_i_i_C(3) + pkin(7) + qJ(4);
	t235 = t253 * t220 - sin(qJ(2)) * pkin(2);
	t244 = t218 * qJD(4);
	t262 = (-t215 * t218 + t235) * qJD(2) + (pkin(4) * sin(pkin(9)) + qJ(3) + pkin(6)) * qJD(1) + t244;
	t221 = pkin(9) + qJ(5);
	t217 = sin(t221);
	t219 = cos(t221);
	t236 = r_i_i_C(1) * t219 - r_i_i_C(2) * t217 + t215;
	t231 = -t236 * t218 + t235;
	t229 = cos(qJ(1));
	t245 = qJD(5) * t220;
	t240 = -qJD(1) + t245;
	t260 = t229 * t240;
	t258 = -t253 * t218 - cos(qJ(2)) * pkin(2);
	t239 = qJD(1) * t220 - qJD(5);
	t227 = sin(qJ(1));
	t248 = qJD(2) * t227;
	t257 = -t218 * t248 + t239 * t229;
	t251 = qJD(1) * t227;
	t250 = qJD(1) * t229;
	t249 = qJD(2) * t220;
	t247 = qJD(2) * t229;
	t246 = qJD(5) * t218;
	t238 = r_i_i_C(1) * t217 + r_i_i_C(2) * t219;
	t237 = t240 * t227;
	t233 = t218 * t247 + t239 * t227;
	t232 = qJD(3) + (-t215 * t220 - pkin(1) + t258) * qJD(1);
	t230 = qJD(4) * t220 + t238 * t246 + (-t236 * t220 + t258) * qJD(2);
	t214 = t217 * t237 - t257 * t219;
	t213 = t257 * t217 + t219 * t237;
	t212 = t217 * t260 + t233 * t219;
	t211 = t233 * t217 - t219 * t260;
	t1 = [t214 * r_i_i_C(1) + t213 * r_i_i_C(2) - t262 * t227 + t232 * t229, t230 * t229 - t231 * t251, t250, -t218 * t251 + t220 * t247, t211 * r_i_i_C(1) + t212 * r_i_i_C(2); -t212 * r_i_i_C(1) + t211 * r_i_i_C(2) + t232 * t227 + t262 * t229, t230 * t227 + t231 * t250, t251, t218 * t250 + t220 * t248, -t213 * r_i_i_C(1) + t214 * r_i_i_C(2); 0, t231 * qJD(2) - t238 * t245 + t244, 0, qJD(2) * t218, (t217 * t246 - t219 * t249) * r_i_i_C(2) + (-t217 * t249 - t219 * t246) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end