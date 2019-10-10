% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPPRR4
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:30
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPPRR4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR4_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPPRR4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPPRR4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:30:55
	% EndTime: 2019-10-09 23:30:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:30:55
	% EndTime: 2019-10-09 23:30:55
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:30:55
	% EndTime: 2019-10-09 23:30:55
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->6), mult. (20->10), div. (0->0), fcn. (12->2), ass. (0->5)
	t10 = -pkin(1) - r_i_i_C(1);
	t9 = r_i_i_C(3) + qJ(2);
	t8 = cos(qJ(1));
	t7 = sin(qJ(1));
	t1 = [t8 * qJD(2) + (t10 * t8 - t9 * t7) * qJD(1), qJD(1) * t8, 0, 0, 0, 0; t7 * qJD(2) + (t10 * t7 + t9 * t8) * qJD(1), qJD(1) * t7, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:30:55
	% EndTime: 2019-10-09 23:30:55
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (14->11), mult. (36->20), div. (0->0), fcn. (26->4), ass. (0->8)
	t58 = -pkin(1) - pkin(2);
	t57 = cos(qJ(1));
	t56 = sin(qJ(1));
	t55 = cos(pkin(9));
	t54 = sin(pkin(9));
	t53 = (t54 * t57 - t55 * t56) * qJD(1);
	t52 = (t54 * t56 + t55 * t57) * qJD(1);
	t1 = [-t52 * r_i_i_C(1) + t53 * r_i_i_C(2) + t57 * qJD(2) + (-qJ(2) * t56 + t57 * t58) * qJD(1), qJD(1) * t57, 0, 0, 0, 0; t53 * r_i_i_C(1) + t52 * r_i_i_C(2) + t56 * qJD(2) + (qJ(2) * t57 + t56 * t58) * qJD(1), qJD(1) * t56, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:30:55
	% EndTime: 2019-10-09 23:30:55
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (28->17), mult. (70->23), div. (0->0), fcn. (58->4), ass. (0->13)
	t31 = -pkin(1) - pkin(2);
	t30 = pkin(3) - r_i_i_C(2);
	t29 = -r_i_i_C(3) - qJ(4);
	t24 = sin(qJ(1));
	t28 = qJD(1) * t24;
	t25 = cos(qJ(1));
	t27 = qJD(1) * t25;
	t22 = sin(pkin(9));
	t23 = cos(pkin(9));
	t26 = t24 * t22 + t25 * t23;
	t20 = t22 * t27 - t23 * t28;
	t19 = t26 * qJD(1);
	t1 = [-t26 * qJD(4) + t25 * qJD(2) + t29 * t20 - t30 * t19 + (-t24 * qJ(2) + t31 * t25) * qJD(1), t27, 0, -t19, 0, 0; -(-t25 * t22 + t24 * t23) * qJD(4) + t24 * qJD(2) + t30 * t20 + t29 * t19 + (t25 * qJ(2) + t31 * t24) * qJD(1), t28, 0, t20, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:30:55
	% EndTime: 2019-10-09 23:30:56
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (59->26), mult. (162->42), div. (0->0), fcn. (140->6), ass. (0->20)
	t61 = -pkin(1) - pkin(2);
	t60 = cos(pkin(9));
	t50 = cos(qJ(1));
	t59 = qJD(1) * t50;
	t47 = sin(qJ(5));
	t58 = qJD(5) * t47;
	t49 = cos(qJ(5));
	t57 = qJD(5) * t49;
	t56 = pkin(3) + pkin(7) + r_i_i_C(3);
	t48 = sin(qJ(1));
	t55 = t48 * t60;
	t54 = -r_i_i_C(1) * t47 - r_i_i_C(2) * t49 - qJ(4);
	t53 = (r_i_i_C(1) * t49 - r_i_i_C(2) * t47) * qJD(5);
	t46 = sin(pkin(9));
	t52 = t48 * t46 + t50 * t60;
	t51 = qJD(4) + t53;
	t42 = t50 * t46 - t55;
	t40 = -qJD(1) * t55 + t46 * t59;
	t39 = t52 * qJD(1);
	t1 = [t50 * qJD(2) - t51 * t52 + t54 * t40 - t56 * t39 + (-qJ(2) * t48 + t61 * t50) * qJD(1), t59, 0, -t39, (t39 * t47 - t42 * t57) * r_i_i_C(2) + (-t39 * t49 - t42 * t58) * r_i_i_C(1), 0; t48 * qJD(2) + t51 * t42 + t56 * t40 + t54 * t39 + (qJ(2) * t50 + t61 * t48) * qJD(1), qJD(1) * t48, 0, t40, (-t40 * t47 - t52 * t57) * r_i_i_C(2) + (t40 * t49 - t52 * t58) * r_i_i_C(1), 0; 0, 0, 0, 0, t53, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:30:57
	% EndTime: 2019-10-09 23:30:57
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (181->49), mult. (506->73), div. (0->0), fcn. (485->8), ass. (0->39)
	t232 = sin(qJ(5));
	t231 = sin(qJ(6));
	t233 = cos(qJ(6));
	t245 = t233 * r_i_i_C(1) - t231 * r_i_i_C(2);
	t241 = pkin(5) + t245;
	t234 = cos(qJ(5));
	t259 = pkin(8) + r_i_i_C(3);
	t248 = t259 * t234;
	t236 = -t241 * t232 + t248;
	t266 = t236 * qJD(5);
	t249 = t259 * t232;
	t237 = t241 * t234 + t249;
	t255 = sin(pkin(9));
	t256 = cos(pkin(9));
	t257 = sin(qJ(1));
	t258 = cos(qJ(1));
	t226 = t258 * t255 - t257 * t256;
	t263 = qJD(1) * t257;
	t262 = qJD(1) * t258;
	t261 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
	t260 = pkin(3) + pkin(7);
	t254 = qJD(5) * t232;
	t253 = qJD(5) * t234;
	t252 = qJD(6) * t226;
	t251 = qJD(6) * t232;
	t250 = qJD(6) * t234;
	t244 = t231 * r_i_i_C(1) + t233 * r_i_i_C(2);
	t225 = -t257 * t255 - t258 * t256;
	t223 = t225 * qJD(1);
	t243 = t225 * t251 - t223;
	t224 = t226 * qJD(1);
	t242 = -t226 * t251 + t224;
	t240 = qJD(6) * t244;
	t239 = -qJD(6) * t225 + t223 * t232 + t226 * t253;
	t238 = -t224 * t232 + t225 * t253 + t252;
	t235 = t237 * qJD(5) - t232 * t240;
	t222 = t242 * t231 + t239 * t233;
	t221 = -t239 * t231 + t242 * t233;
	t1 = [t245 * t252 - (-t244 - t260) * t223 - qJ(2) * t263 + (-qJ(4) + t236) * t224 + (qJD(4) + t235) * t225 + t261 * t258, t262, 0, t223, t237 * t223 + (-t234 * t240 + t266) * t226, t221 * r_i_i_C(1) - t222 * r_i_i_C(2); t222 * r_i_i_C(1) + t221 * r_i_i_C(2) + t260 * t224 + qJ(2) * t262 + (qJD(4) + (pkin(5) * t234 + t249) * qJD(5)) * t226 - (-pkin(5) * t232 - qJ(4) + t248) * t223 + t261 * t257, t263, 0, t224, t237 * t224 + (t244 * t250 - t266) * t225, (t243 * r_i_i_C(1) + t238 * r_i_i_C(2)) * t233 + (t238 * r_i_i_C(1) - t243 * r_i_i_C(2)) * t231; 0, 0, 0, 0, t235, (-t231 * t250 - t233 * t254) * r_i_i_C(2) + (-t231 * t254 + t233 * t250) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end