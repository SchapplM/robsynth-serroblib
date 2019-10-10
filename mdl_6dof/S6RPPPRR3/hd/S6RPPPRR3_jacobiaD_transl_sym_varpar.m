% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPPRR3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:29
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPPRR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR3_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPPRR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPPRR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
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
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
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
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (14->11), mult. (36->20), div. (0->0), fcn. (26->4), ass. (0->8)
	t58 = -pkin(1) - pkin(2);
	t57 = cos(qJ(1));
	t56 = sin(qJ(1));
	t55 = cos(pkin(9));
	t54 = sin(pkin(9));
	t53 = (t54 * t57 - t55 * t56) * qJD(1);
	t52 = (t54 * t56 + t55 * t57) * qJD(1);
	t1 = [-t52 * r_i_i_C(1) + t53 * r_i_i_C(2) + t57 * qJD(2) + (-qJ(2) * t56 + t58 * t57) * qJD(1), qJD(1) * t57, 0, 0, 0, 0; t53 * r_i_i_C(1) + t52 * r_i_i_C(2) + t56 * qJD(2) + (qJ(2) * t57 + t58 * t56) * qJD(1), qJD(1) * t56, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (32->17), mult. (84->24), div. (0->0), fcn. (70->6), ass. (0->12)
	t36 = -pkin(1) - pkin(2);
	t35 = -r_i_i_C(3) - qJ(4);
	t27 = sin(pkin(9));
	t29 = cos(pkin(9));
	t30 = sin(qJ(1));
	t31 = cos(qJ(1));
	t34 = t31 * t27 - t30 * t29;
	t33 = t30 * t27 + t31 * t29;
	t32 = r_i_i_C(1) * cos(pkin(10)) - r_i_i_C(2) * sin(pkin(10)) + pkin(3);
	t25 = t34 * qJD(1);
	t24 = t33 * qJD(1);
	t1 = [-t33 * qJD(4) + t31 * qJD(2) + t35 * t25 - t32 * t24 + (-t30 * qJ(2) + t36 * t31) * qJD(1), qJD(1) * t31, 0, -t24, 0, 0; t34 * qJD(4) + t30 * qJD(2) + t35 * t24 + t32 * t25 + (t31 * qJ(2) + t36 * t30) * qJD(1), qJD(1) * t30, 0, t25, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (77->28), mult. (154->44), div. (0->0), fcn. (134->7), ass. (0->18)
	t52 = pkin(10) + qJ(5);
	t50 = sin(t52);
	t51 = cos(t52);
	t58 = (r_i_i_C(1) * t50 + r_i_i_C(2) * t51) * qJD(5);
	t64 = -pkin(1) - pkin(2);
	t63 = -r_i_i_C(3) - pkin(7) - qJ(4);
	t62 = qJD(5) * t50;
	t61 = qJD(5) * t51;
	t53 = sin(pkin(9));
	t54 = cos(pkin(9));
	t56 = sin(qJ(1));
	t57 = cos(qJ(1));
	t46 = t53 * t57 - t56 * t54;
	t47 = t56 * t53 + t54 * t57;
	t59 = r_i_i_C(1) * t51 - r_i_i_C(2) * t50 + cos(pkin(10)) * pkin(4) + pkin(3);
	t45 = t46 * qJD(1);
	t44 = t47 * qJD(1);
	t1 = [qJD(2) * t57 - t47 * qJD(4) + t63 * t45 - t46 * t58 - t59 * t44 + (-qJ(2) * t56 + t64 * t57) * qJD(1), qJD(1) * t57, 0, -t44, (-t45 * t51 + t47 * t62) * r_i_i_C(2) + (-t45 * t50 - t47 * t61) * r_i_i_C(1), 0; t56 * qJD(2) + t46 * qJD(4) + t63 * t44 - t47 * t58 + t59 * t45 + (qJ(2) * t57 + t64 * t56) * qJD(1), qJD(1) * t56, 0, t45, (-t44 * t51 - t46 * t62) * r_i_i_C(2) + (-t44 * t50 + t46 * t61) * r_i_i_C(1), 0; 0, 0, 0, 0, t58, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:16
	% EndTime: 2019-10-09 23:29:17
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (252->53), mult. (498->82), div. (0->0), fcn. (479->9), ass. (0->40)
	t245 = pkin(10) + qJ(5);
	t243 = sin(t245);
	t247 = sin(qJ(6));
	t248 = cos(qJ(6));
	t254 = t248 * r_i_i_C(1) - t247 * r_i_i_C(2) + pkin(5);
	t244 = cos(t245);
	t271 = pkin(8) + r_i_i_C(3);
	t261 = t271 * t244;
	t250 = -t254 * t243 + t261;
	t278 = qJD(5) * t250;
	t260 = t271 * t243;
	t277 = t254 * t244 + t260;
	t267 = sin(pkin(9));
	t268 = cos(pkin(9));
	t269 = sin(qJ(1));
	t270 = cos(qJ(1));
	t237 = t270 * t267 - t269 * t268;
	t274 = qJD(1) * t269;
	t273 = qJD(1) * t270;
	t272 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
	t266 = qJD(5) * t243;
	t265 = qJD(5) * t244;
	t264 = qJD(6) * t244;
	t263 = qJD(6) * t247;
	t262 = qJD(6) * t248;
	t257 = t247 * r_i_i_C(1) + t248 * r_i_i_C(2);
	t236 = -t267 * t269 - t268 * t270;
	t234 = t236 * qJD(1);
	t256 = t236 * t264 + t234;
	t235 = t237 * qJD(1);
	t255 = t237 * t264 + t235;
	t253 = qJD(6) * t257;
	t252 = qJD(6) * t236 + t234 * t244 - t237 * t266;
	t251 = qJD(6) * t237 + t235 * t244 + t236 * t266;
	t249 = t277 * qJD(5) - t243 * t253;
	t246 = -pkin(7) - qJ(4);
	t242 = cos(pkin(10)) * pkin(4) + pkin(3);
	t233 = t247 * t256 + t248 * t251;
	t232 = -t247 * t251 + t248 * t256;
	t1 = [(-t235 * t247 + t236 * t262) * r_i_i_C(1) + (-t235 * t248 - t236 * t263) * r_i_i_C(2) + t235 * t246 + t236 * qJD(4) - qJ(2) * t274 - (-t242 - t277) * t234 + (-t244 * t253 + t278) * t237 + t272 * t270, t273, 0, t234, t235 * t250 + t236 * t249, t232 * r_i_i_C(1) - t233 * r_i_i_C(2); t233 * r_i_i_C(1) + t232 * r_i_i_C(2) + t237 * qJD(4) - t234 * t246 + (pkin(5) * t243 - t261) * t236 * qJD(5) + qJ(2) * t273 + (pkin(5) * t244 + t242 + t260) * t235 + t272 * t269, t274, 0, t235, -t234 * t250 + t237 * t249, (r_i_i_C(1) * t255 + r_i_i_C(2) * t252) * t248 + (r_i_i_C(1) * t252 - r_i_i_C(2) * t255) * t247; 0, 0, 0, 0, t257 * t264 - t278, (-t243 * t263 + t248 * t265) * r_i_i_C(2) + (t243 * t262 + t247 * t265) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end