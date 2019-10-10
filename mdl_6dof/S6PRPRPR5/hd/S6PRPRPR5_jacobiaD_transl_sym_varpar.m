% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPRPR5
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:37
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRPR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR5_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->7), mult. (30->20), div. (0->0), fcn. (24->6), ass. (0->8)
	t53 = cos(pkin(6));
	t54 = sin(qJ(2));
	t57 = t53 * t54;
	t55 = cos(qJ(2));
	t56 = t53 * t55;
	t52 = cos(pkin(10));
	t50 = sin(pkin(10));
	t1 = [0, ((t50 * t57 - t52 * t55) * r_i_i_C(1) + (t50 * t56 + t52 * t54) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, ((-t50 * t55 - t52 * t57) * r_i_i_C(1) + (t50 * t54 - t52 * t56) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t54 - r_i_i_C(2) * t55) * sin(pkin(6)) * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:22
	% EndTime: 2019-10-09 21:37:22
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (29->14), mult. (103->28), div. (0->0), fcn. (90->8), ass. (0->16)
	t153 = sin(pkin(10));
	t156 = cos(pkin(10));
	t159 = cos(qJ(2));
	t157 = cos(pkin(6));
	t158 = sin(qJ(2));
	t165 = t157 * t158;
	t168 = t153 * t165 - t156 * t159;
	t167 = r_i_i_C(3) + qJ(3);
	t164 = t157 * t159;
	t162 = qJD(2) * t167;
	t161 = -cos(pkin(11)) * r_i_i_C(1) + sin(pkin(11)) * r_i_i_C(2) - pkin(2);
	t160 = t153 * t159 + t156 * t165;
	t154 = sin(pkin(6));
	t150 = t168 * qJD(2);
	t148 = t160 * qJD(2);
	t1 = [0, -t168 * qJD(3) - (t153 * t164 + t156 * t158) * t162 - t161 * t150, -t150, 0, 0, 0; 0, t160 * qJD(3) - (t153 * t158 - t156 * t164) * t162 + t161 * t148, t148, 0, 0, 0; 0, (qJD(3) * t158 + (t161 * t158 + t167 * t159) * qJD(2)) * t154, t154 * qJD(2) * t158, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:22
	% EndTime: 2019-10-09 21:37:22
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (98->33), mult. (216->67), div. (0->0), fcn. (199->9), ass. (0->30)
	t212 = r_i_i_C(3) + pkin(8) + qJ(3);
	t192 = sin(pkin(10));
	t193 = sin(pkin(6));
	t211 = t192 * t193;
	t194 = cos(pkin(10));
	t210 = t193 * t194;
	t197 = sin(qJ(2));
	t209 = t193 * t197;
	t195 = cos(pkin(6));
	t208 = t195 * t197;
	t198 = cos(qJ(2));
	t207 = t195 * t198;
	t206 = qJD(2) * t197;
	t205 = qJD(2) * t198;
	t204 = t192 * t206;
	t203 = t194 * t205;
	t191 = pkin(11) + qJ(4);
	t189 = sin(t191);
	t190 = cos(t191);
	t202 = t189 * r_i_i_C(1) + t190 * r_i_i_C(2);
	t201 = -t190 * r_i_i_C(1) + t189 * r_i_i_C(2) - cos(pkin(11)) * pkin(3) - pkin(2);
	t183 = t192 * t198 + t194 * t208;
	t200 = t192 * t207 + t194 * t197;
	t199 = qJD(4) * t202;
	t185 = -t192 * t208 + t194 * t198;
	t181 = -t195 * t204 + t203;
	t180 = t200 * qJD(2);
	t179 = t183 * qJD(2);
	t178 = -t195 * t203 + t204;
	t1 = [0, t185 * qJD(3) - t212 * t180 + t201 * t181 + t200 * t199, t181, t202 * t180 + ((-t185 * t190 - t189 * t211) * r_i_i_C(1) + (t185 * t189 - t190 * t211) * r_i_i_C(2)) * qJD(4), 0, 0; 0, t183 * qJD(3) - t212 * t178 - (-t192 * t197 + t194 * t207) * t199 + t201 * t179, t179, t202 * t178 + ((-t183 * t190 + t189 * t210) * r_i_i_C(1) + (t183 * t189 + t190 * t210) * r_i_i_C(2)) * qJD(4), 0, 0; 0, (qJD(3) * t197 - t198 * t199 + (t201 * t197 + t212 * t198) * qJD(2)) * t193, t193 * t206, -t202 * t193 * t205 + ((-t189 * t195 - t190 * t209) * r_i_i_C(1) + (t189 * t209 - t190 * t195) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:23
	% EndTime: 2019-10-09 21:37:23
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (215->42), mult. (429->77), div. (0->0), fcn. (411->9), ass. (0->38)
	t261 = -pkin(4) + r_i_i_C(2);
	t260 = r_i_i_C(1) + pkin(8) + qJ(3);
	t259 = r_i_i_C(3) + qJ(5);
	t236 = sin(pkin(10));
	t237 = sin(pkin(6));
	t258 = t236 * t237;
	t238 = cos(pkin(10));
	t257 = t237 * t238;
	t241 = sin(qJ(2));
	t256 = t237 * t241;
	t239 = cos(pkin(6));
	t255 = t239 * t241;
	t242 = cos(qJ(2));
	t254 = t239 * t242;
	t253 = qJD(2) * t241;
	t252 = qJD(2) * t242;
	t251 = t237 * t252;
	t250 = t236 * t253;
	t249 = t238 * t252;
	t227 = t236 * t242 + t238 * t255;
	t235 = pkin(11) + qJ(4);
	t233 = sin(t235);
	t234 = cos(t235);
	t248 = -t227 * t234 + t233 * t257;
	t229 = -t236 * t255 + t238 * t242;
	t247 = t229 * t234 + t233 * t258;
	t246 = t239 * t233 + t234 * t256;
	t245 = t236 * t254 + t238 * t241;
	t244 = -t259 * t233 + t261 * t234 - cos(pkin(11)) * pkin(3) - pkin(2);
	t243 = qJD(5) * t233 + (t261 * t233 + t259 * t234) * qJD(4);
	t225 = -t239 * t250 + t249;
	t224 = t245 * qJD(2);
	t223 = t227 * qJD(2);
	t222 = -t239 * t249 + t250;
	t220 = t246 * qJD(4) + t233 * t251;
	t218 = t247 * qJD(4) - t224 * t233;
	t216 = -t248 * qJD(4) - t222 * t233;
	t1 = [0, t229 * qJD(3) - t260 * t224 + t244 * t225 - t243 * t245, t225, t247 * qJD(5) + t259 * (-t224 * t234 + (-t229 * t233 + t234 * t258) * qJD(4)) + t261 * t218, t218, 0; 0, t227 * qJD(3) - t260 * t222 + t243 * (-t236 * t241 + t238 * t254) + t244 * t223, t223, -t248 * qJD(5) + t259 * (-t222 * t234 + (-t227 * t233 - t234 * t257) * qJD(4)) + t261 * t216, t216, 0; 0, (qJD(3) * t241 + t243 * t242 + (t244 * t241 + t260 * t242) * qJD(2)) * t237, t237 * t253, t246 * qJD(5) + t259 * (t234 * t251 + (-t233 * t256 + t234 * t239) * qJD(4)) + t261 * t220, t220, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:23
	% EndTime: 2019-10-09 21:37:24
	% DurationCPUTime: 0.51s
	% Computational Cost: add. (444->66), mult. (905->120), div. (0->0), fcn. (907->11), ass. (0->48)
	t311 = pkin(11) + qJ(4);
	t309 = sin(t311);
	t310 = cos(t311);
	t315 = sin(qJ(6));
	t317 = cos(qJ(6));
	t332 = t317 * r_i_i_C(1) - t315 * r_i_i_C(2);
	t321 = t332 * qJD(6) + qJD(5);
	t331 = -t315 * r_i_i_C(1) - t317 * r_i_i_C(2);
	t329 = qJ(5) - t331;
	t340 = pkin(4) + pkin(9) + r_i_i_C(3);
	t319 = (t340 * t309 - t329 * t310) * qJD(4) - t321 * t309;
	t312 = sin(pkin(10));
	t316 = sin(qJ(2));
	t318 = cos(qJ(2));
	t346 = cos(pkin(10));
	t347 = cos(pkin(6));
	t330 = t347 * t346;
	t301 = t312 * t318 + t316 * t330;
	t313 = sin(pkin(6));
	t336 = t313 * t346;
	t349 = t301 * t310 - t309 * t336;
	t337 = t312 * t347;
	t303 = -t316 * t337 + t346 * t318;
	t344 = t312 * t313;
	t343 = t313 * t316;
	t342 = t313 * t318;
	t341 = qJD(2) * t316;
	t339 = qJD(2) * t342;
	t338 = t313 * t341;
	t328 = t318 * t330;
	t327 = -t303 * t309 + t310 * t344;
	t326 = t303 * t310 + t309 * t344;
	t325 = pkin(5) + pkin(8) + qJ(3) + t332;
	t294 = t309 * t343 - t347 * t310;
	t324 = t347 * t309 + t310 * t343;
	t323 = -t301 * t309 - t310 * t336;
	t322 = t331 * qJD(6) + qJD(3);
	t302 = t346 * t316 + t318 * t337;
	t320 = -t329 * t309 - t340 * t310 - cos(pkin(11)) * pkin(3) - pkin(2);
	t300 = t312 * t316 - t328;
	t299 = t303 * qJD(2);
	t298 = t302 * qJD(2);
	t297 = t301 * qJD(2);
	t296 = -qJD(2) * t328 + t312 * t341;
	t288 = t324 * qJD(4) + t309 * t339;
	t286 = t326 * qJD(4) - t298 * t309;
	t284 = t349 * qJD(4) - t296 * t309;
	t1 = [0, -t325 * t298 + t320 * t299 + t319 * t302 + t322 * t303, t299, t321 * t326 + t329 * (t327 * qJD(4) - t298 * t310) - t340 * t286, t286, (t286 * t317 - t299 * t315) * r_i_i_C(1) + (-t286 * t315 - t299 * t317) * r_i_i_C(2) + ((-t302 * t317 + t315 * t327) * r_i_i_C(1) + (t302 * t315 + t317 * t327) * r_i_i_C(2)) * qJD(6); 0, -t325 * t296 + t320 * t297 + t319 * t300 + t322 * t301, t297, t321 * t349 + t329 * (t323 * qJD(4) - t296 * t310) - t340 * t284, t284, (t284 * t317 - t297 * t315) * r_i_i_C(1) + (-t284 * t315 - t297 * t317) * r_i_i_C(2) + ((-t300 * t317 + t315 * t323) * r_i_i_C(1) + (t300 * t315 + t317 * t323) * r_i_i_C(2)) * qJD(6); 0, ((t320 * qJD(2) + t322) * t316 + (t325 * qJD(2) - t319) * t318) * t313, t338, t321 * t324 + t329 * (-t294 * qJD(4) + t310 * t339) - t340 * t288, t288, (t288 * t317 - t315 * t338) * r_i_i_C(1) + (-t288 * t315 - t317 * t338) * r_i_i_C(2) + ((-t294 * t315 + t317 * t342) * r_i_i_C(1) + (-t294 * t317 - t315 * t342) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end