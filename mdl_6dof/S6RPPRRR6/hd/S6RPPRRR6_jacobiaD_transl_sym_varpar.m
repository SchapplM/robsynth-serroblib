% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:10
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRR6_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR6_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR6_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR6_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:10:02
	% EndTime: 2019-10-10 00:10:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:10:02
	% EndTime: 2019-10-10 00:10:02
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
	% StartTime: 2019-10-10 00:10:02
	% EndTime: 2019-10-10 00:10:02
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->6), mult. (20->10), div. (0->0), fcn. (12->2), ass. (0->5)
	t10 = -pkin(1) + r_i_i_C(2);
	t9 = r_i_i_C(3) + qJ(2);
	t8 = cos(qJ(1));
	t7 = sin(qJ(1));
	t1 = [t8 * qJD(2) + (t10 * t8 - t9 * t7) * qJD(1), qJD(1) * t8, 0, 0, 0, 0; t7 * qJD(2) + (t10 * t7 + t9 * t8) * qJD(1), qJD(1) * t7, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:10:02
	% EndTime: 2019-10-10 00:10:02
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (13->9), mult. (28->12), div. (0->0), fcn. (18->2), ass. (0->7)
	t12 = r_i_i_C(2) + qJ(2);
	t8 = sin(qJ(1));
	t11 = qJD(1) * t8;
	t10 = -pkin(1) - r_i_i_C(3) - qJ(3);
	t9 = cos(qJ(1));
	t7 = qJD(1) * t9;
	t1 = [t9 * qJD(2) - t8 * qJD(3) + (t10 * t9 - t12 * t8) * qJD(1), t7, -t11, 0, 0, 0; t8 * qJD(2) + t9 * qJD(3) + (t10 * t8 + t12 * t9) * qJD(1), t11, t7, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:10:02
	% EndTime: 2019-10-10 00:10:02
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (28->20), mult. (80->31), div. (0->0), fcn. (52->4), ass. (0->13)
	t18 = sin(qJ(4));
	t20 = cos(qJ(4));
	t23 = (r_i_i_C(1) * t20 - r_i_i_C(2) * t18) * qJD(4);
	t29 = qJD(3) + t23;
	t19 = sin(qJ(1));
	t28 = qJD(1) * t19;
	t21 = cos(qJ(1));
	t17 = qJD(1) * t21;
	t27 = qJD(4) * t19;
	t26 = qJD(4) * t21;
	t25 = pkin(7) + r_i_i_C(3) - qJ(2);
	t22 = -r_i_i_C(1) * t18 - r_i_i_C(2) * t20 - pkin(1) - qJ(3);
	t1 = [t21 * qJD(2) - t29 * t19 + (t19 * t25 + t21 * t22) * qJD(1), t17, -t28, (t18 * t28 - t20 * t26) * r_i_i_C(2) + (-t18 * t26 - t20 * t28) * r_i_i_C(1), 0, 0; t19 * qJD(2) + t29 * t21 + (t19 * t22 - t21 * t25) * qJD(1), t28, t17, (-t17 * t18 - t20 * t27) * r_i_i_C(2) + (t17 * t20 - t18 * t27) * r_i_i_C(1), 0, 0; 0, 0, 0, -t23, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:10:03
	% EndTime: 2019-10-10 00:10:03
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (92->44), mult. (286->73), div. (0->0), fcn. (223->6), ass. (0->33)
	t206 = cos(qJ(4));
	t203 = sin(qJ(4));
	t228 = pkin(8) + r_i_i_C(3);
	t229 = t228 * t203;
	t233 = (pkin(4) * t206 + t229) * qJD(4) + qJD(3);
	t202 = sin(qJ(5));
	t205 = cos(qJ(5));
	t212 = r_i_i_C(1) * t205 - r_i_i_C(2) * t202 + pkin(4);
	t231 = t212 * t206 + t229;
	t226 = pkin(7) - qJ(2);
	t207 = cos(qJ(1));
	t225 = t205 * t207;
	t204 = sin(qJ(1));
	t224 = qJD(1) * t204;
	t201 = qJD(1) * t207;
	t223 = qJD(4) * t203;
	t222 = qJD(4) * t206;
	t221 = qJD(4) * t207;
	t220 = qJD(5) * t203;
	t219 = qJD(5) * t206;
	t216 = t228 * t206;
	t215 = qJD(1) + t220;
	t214 = qJD(1) * t203 + qJD(5);
	t213 = r_i_i_C(1) * t202 + r_i_i_C(2) * t205;
	t211 = t215 * t202;
	t210 = t213 * qJD(5);
	t209 = -pkin(4) * t203 - pkin(1) - qJ(3) + t216;
	t208 = t214 * t204 - t206 * t221;
	t200 = -t214 * t225 + (-t205 * t222 + t211) * t204;
	t199 = t215 * t205 * t204 + (t204 * t222 + t214 * t207) * t202;
	t198 = t208 * t205 + t207 * t211;
	t197 = t208 * t202 - t215 * t225;
	t1 = [t200 * r_i_i_C(1) + t199 * r_i_i_C(2) + t207 * qJD(2) - t233 * t204 + (t226 * t204 + t209 * t207) * qJD(1), t201, -t224, (-t212 * t221 - t228 * t224) * t203 + (-t212 * t224 + (t228 * qJD(4) - t210) * t207) * t206, t197 * r_i_i_C(1) + t198 * r_i_i_C(2), 0; -t198 * r_i_i_C(1) + t197 * r_i_i_C(2) + t204 * qJD(2) + t233 * t207 + (t209 * t204 - t226 * t207) * qJD(1), t224, t201, t231 * t201 + (-t206 * t210 + (-t212 * t203 + t216) * qJD(4)) * t204, -t199 * r_i_i_C(1) + t200 * r_i_i_C(2), 0; 0, 0, 0, -t231 * qJD(4) + t213 * t220, (t202 * t219 + t205 * t223) * r_i_i_C(2) + (t202 * t223 - t205 * t219) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:10:03
	% EndTime: 2019-10-10 00:10:03
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (270->57), mult. (436->85), div. (0->0), fcn. (344->8), ass. (0->49)
	t240 = cos(qJ(5));
	t231 = t240 * pkin(5) + pkin(4);
	t238 = sin(qJ(4));
	t241 = cos(qJ(4));
	t237 = sin(qJ(5));
	t278 = pkin(5) * t237;
	t263 = qJD(5) * t278;
	t273 = r_i_i_C(3) + pkin(9) + pkin(8);
	t281 = t273 * t238;
	t288 = (t231 * t241 + t281) * qJD(4) - t238 * t263 + qJD(3);
	t236 = qJ(5) + qJ(6);
	t233 = sin(t236);
	t234 = cos(t236);
	t287 = r_i_i_C(1) * t233 + r_i_i_C(2) * t234;
	t286 = -r_i_i_C(1) * t234 + r_i_i_C(2) * t233;
	t248 = t231 - t286;
	t284 = t248 * t241 + t281;
	t242 = cos(qJ(1));
	t235 = qJD(5) + qJD(6);
	t257 = t235 * t238 + qJD(1);
	t283 = t242 * t257;
	t280 = t287 * t235 + t263;
	t268 = qJD(1) * t238;
	t256 = t235 + t268;
	t239 = sin(qJ(1));
	t262 = qJD(4) * t239 * t241;
	t279 = t256 * t242 + t262;
	t265 = qJD(4) * t242;
	t261 = t241 * t265;
	t245 = t256 * t239 - t261;
	t224 = t245 * t233 - t234 * t283;
	t225 = t233 * t283 + t245 * t234;
	t270 = t224 * r_i_i_C(1) + t225 * r_i_i_C(2);
	t251 = t257 * t239;
	t226 = t279 * t233 + t234 * t251;
	t227 = t233 * t251 - t279 * t234;
	t269 = -t226 * r_i_i_C(1) + t227 * r_i_i_C(2);
	t267 = qJD(1) * t239;
	t232 = qJD(1) * t242;
	t266 = qJD(4) * t238;
	t264 = qJD(5) * t240;
	t258 = t273 * t241;
	t255 = qJD(5) + t268;
	t254 = pkin(7) - qJ(2) + t278;
	t253 = -pkin(5) * t264 + qJD(2);
	t250 = (-qJD(5) * t238 - qJD(1)) * t240;
	t247 = t286 * t235 * t241 + t287 * t266;
	t246 = -t231 * t238 - pkin(1) - qJ(3) + t258;
	t1 = [t227 * r_i_i_C(1) + t226 * r_i_i_C(2) + t253 * t242 - t288 * t239 + (t254 * t239 + t246 * t242) * qJD(1), t232, -t267, (-t248 * t265 - t273 * t267) * t238 + (-t248 * t267 + (t273 * qJD(4) - t280) * t242) * t241, (t242 * t250 + (t255 * t239 - t261) * t237) * pkin(5) + t270, t270; -t225 * r_i_i_C(1) + t224 * r_i_i_C(2) + t253 * t239 + t288 * t242 + (t246 * t239 - t254 * t242) * qJD(1), t267, t232, t284 * t232 + (-t280 * t241 + (-t248 * t238 + t258) * qJD(4)) * t239, (t239 * t250 + (-t255 * t242 - t262) * t237) * pkin(5) + t269, t269; 0, 0, 0, -t284 * qJD(4) + t280 * t238, (t237 * t266 - t241 * t264) * pkin(5) + t247, t247;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end