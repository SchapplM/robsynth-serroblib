% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:10
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPPR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR3_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPPR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
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
	% StartTime: 2019-10-09 22:10:55
	% EndTime: 2019-10-09 22:10:55
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (55->26), mult. (191->57), div. (0->0), fcn. (172->8), ass. (0->24)
	t180 = sin(pkin(10));
	t182 = cos(pkin(10));
	t185 = sin(qJ(2));
	t183 = cos(pkin(6));
	t187 = cos(qJ(2));
	t195 = t183 * t187;
	t201 = -t180 * t185 + t182 * t195;
	t200 = pkin(8) + r_i_i_C(3);
	t181 = sin(pkin(6));
	t184 = sin(qJ(3));
	t198 = t181 * t184;
	t186 = cos(qJ(3));
	t197 = t181 * t186;
	t196 = t183 * t185;
	t193 = t184 * r_i_i_C(1) + t186 * r_i_i_C(2);
	t192 = t186 * r_i_i_C(1) - t184 * r_i_i_C(2) + pkin(2);
	t176 = t180 * t187 + t182 * t196;
	t191 = t180 * t195 + t182 * t185;
	t190 = t180 * t196 - t182 * t187;
	t189 = qJD(3) * t193;
	t188 = qJD(2) * t192;
	t173 = t191 * qJD(2);
	t171 = t201 * qJD(2);
	t1 = [0, -t200 * t173 + t190 * t188 + t191 * t189, t193 * t173 + ((-t180 * t198 + t186 * t190) * r_i_i_C(1) + (-t180 * t197 - t184 * t190) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, t200 * t171 - t176 * t188 - t201 * t189, -t193 * t171 + ((-t176 * t186 + t182 * t198) * r_i_i_C(1) + (t176 * t184 + t182 * t197) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, (-t187 * t189 + (-t192 * t185 + t200 * t187) * qJD(2)) * t181, -t193 * t187 * t181 * qJD(2) + ((-t183 * t184 - t185 * t197) * r_i_i_C(1) + (-t183 * t186 + t185 * t198) * r_i_i_C(2)) * qJD(3), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:56
	% EndTime: 2019-10-09 22:10:56
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (124->35), mult. (404->66), div. (0->0), fcn. (384->8), ass. (0->32)
	t231 = sin(pkin(10));
	t233 = cos(pkin(10));
	t236 = sin(qJ(2));
	t234 = cos(pkin(6));
	t238 = cos(qJ(2));
	t248 = t234 * t238;
	t259 = -t231 * t236 + t233 * t248;
	t249 = t234 * t236;
	t226 = t231 * t238 + t233 * t249;
	t237 = cos(qJ(3));
	t232 = sin(pkin(6));
	t235 = sin(qJ(3));
	t251 = t232 * t235;
	t258 = -t226 * t237 + t233 * t251;
	t254 = r_i_i_C(3) + qJ(4);
	t256 = pkin(3) + r_i_i_C(1);
	t257 = t254 * t235 + t256 * t237 + pkin(2);
	t255 = pkin(8) + r_i_i_C(2);
	t250 = t232 * t237;
	t245 = qJD(2) * t232 * t238;
	t242 = t231 * t249 - t233 * t238;
	t244 = t231 * t251 - t237 * t242;
	t243 = t231 * t248 + t233 * t236;
	t241 = t234 * t235 + t236 * t250;
	t240 = qJD(2) * t257;
	t239 = qJD(4) * t235 + (-t256 * t235 + t254 * t237) * qJD(3);
	t223 = t243 * qJD(2);
	t221 = t259 * qJD(2);
	t219 = t241 * qJD(3) + t235 * t245;
	t217 = t244 * qJD(3) - t223 * t235;
	t215 = -t258 * qJD(3) + t221 * t235;
	t1 = [0, -t255 * t223 - t239 * t243 + t242 * t240, t244 * qJD(4) + t254 * (-t223 * t237 + (t231 * t250 + t235 * t242) * qJD(3)) - t256 * t217, t217, 0, 0; 0, t255 * t221 - t226 * t240 + t239 * t259, -t258 * qJD(4) + t254 * (t221 * t237 + (-t226 * t235 - t233 * t250) * qJD(3)) - t256 * t215, t215, 0, 0; 0, (t239 * t238 + (-t257 * t236 + t255 * t238) * qJD(2)) * t232, t241 * qJD(4) + t254 * (t237 * t245 + (t234 * t237 - t236 * t251) * qJD(3)) - t256 * t219, t219, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:55
	% EndTime: 2019-10-09 22:10:55
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (160->41), mult. (514->73), div. (0->0), fcn. (488->8), ass. (0->34)
	t197 = sin(pkin(10));
	t199 = cos(pkin(10));
	t204 = cos(qJ(2));
	t200 = cos(pkin(6));
	t202 = sin(qJ(2));
	t218 = t200 * t202;
	t192 = t197 * t204 + t199 * t218;
	t203 = cos(qJ(3));
	t198 = sin(pkin(6));
	t201 = sin(qJ(3));
	t220 = t198 * t201;
	t223 = -t192 * t203 + t199 * t220;
	t215 = pkin(3) + pkin(4) - r_i_i_C(2);
	t222 = r_i_i_C(1) + qJ(4);
	t206 = t222 * t201 + t215 * t203 + pkin(2);
	t219 = t198 * t203;
	t217 = t200 * t204;
	t216 = qJD(2) * t202;
	t214 = pkin(8) - r_i_i_C(3) - qJ(5);
	t212 = t199 * t217;
	t211 = qJD(2) * t198 * t204;
	t208 = t197 * t218 - t199 * t204;
	t210 = t197 * t220 - t203 * t208;
	t209 = t197 * t217 + t199 * t202;
	t207 = t200 * t201 + t202 * t219;
	t205 = qJD(4) * t201 + (-t215 * t201 + t222 * t203) * qJD(3);
	t190 = t208 * qJD(2);
	t189 = t209 * qJD(2);
	t188 = t192 * qJD(2);
	t187 = -qJD(2) * t212 + t197 * t216;
	t185 = t207 * qJD(3) + t201 * t211;
	t183 = t210 * qJD(3) - t189 * t201;
	t181 = -t223 * qJD(3) - t187 * t201;
	t1 = [0, qJD(5) * t208 - t214 * t189 + t206 * t190 - t205 * t209, t210 * qJD(4) + t222 * (-t189 * t203 + (t197 * t219 + t201 * t208) * qJD(3)) - t215 * t183, t183, t190, 0; 0, -t192 * qJD(5) - t214 * t187 + t205 * (-t197 * t202 + t212) - t206 * t188, -t223 * qJD(4) + t222 * (-t187 * t203 + (-t192 * t201 - t199 * t219) * qJD(3)) - t215 * t181, t181, -t188, 0; 0, (-qJD(5) * t202 + t205 * t204 + (-t206 * t202 + t214 * t204) * qJD(2)) * t198, t207 * qJD(4) + t222 * (t203 * t211 + (t200 * t203 - t202 * t220) * qJD(3)) - t215 * t185, t185, -t198 * t216, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:56
	% EndTime: 2019-10-09 22:10:56
	% DurationCPUTime: 0.46s
	% Computational Cost: add. (329->66), mult. (1048->118), div. (0->0), fcn. (1040->10), ass. (0->46)
	t296 = sin(pkin(10));
	t298 = cos(pkin(10));
	t305 = cos(qJ(2));
	t299 = cos(pkin(6));
	t302 = sin(qJ(2));
	t326 = t299 * t302;
	t288 = t296 * t305 + t298 * t326;
	t304 = cos(qJ(3));
	t297 = sin(pkin(6));
	t301 = sin(qJ(3));
	t329 = t297 * t301;
	t331 = t288 * t304 - t298 * t329;
	t300 = sin(qJ(6));
	t303 = cos(qJ(6));
	t318 = -t303 * r_i_i_C(1) + t300 * r_i_i_C(2);
	t311 = pkin(5) + qJ(4) - t318;
	t323 = pkin(3) + pkin(4) + pkin(9) + r_i_i_C(3);
	t307 = t311 * t301 + t323 * t304 + pkin(2);
	t328 = t297 * t304;
	t327 = t297 * t305;
	t325 = t299 * t305;
	t324 = qJD(2) * t302;
	t321 = t298 * t325;
	t320 = qJD(2) * t327;
	t319 = t297 * t324;
	t317 = -t300 * r_i_i_C(1) - t303 * r_i_i_C(2);
	t277 = t288 * t301 + t298 * t328;
	t313 = t296 * t326 - t298 * t305;
	t316 = t296 * t328 + t301 * t313;
	t315 = t296 * t329 - t304 * t313;
	t314 = t296 * t325 + t298 * t302;
	t312 = t299 * t301 + t302 * t328;
	t291 = -t299 * t304 + t302 * t329;
	t310 = -pkin(8) + qJ(5) - t317;
	t309 = t317 * qJD(6) + qJD(4);
	t308 = t318 * qJD(6) - qJD(5);
	t306 = t309 * t301 + (-t323 * t301 + t311 * t304) * qJD(3);
	t287 = -t296 * t302 + t321;
	t286 = t313 * qJD(2);
	t285 = t314 * qJD(2);
	t284 = t288 * qJD(2);
	t283 = -qJD(2) * t321 + t296 * t324;
	t281 = t312 * qJD(3) + t301 * t320;
	t275 = t315 * qJD(3) - t285 * t301;
	t273 = t331 * qJD(3) - t283 * t301;
	t1 = [0, t310 * t285 + t307 * t286 - t306 * t314 - t308 * t313, t309 * t315 + t311 * (t316 * qJD(3) - t285 * t304) - t323 * t275, t275, t286, (-t275 * t300 + t286 * t303) * r_i_i_C(1) + (-t275 * t303 - t286 * t300) * r_i_i_C(2) + ((t300 * t314 + t303 * t316) * r_i_i_C(1) + (-t300 * t316 + t303 * t314) * r_i_i_C(2)) * qJD(6); 0, t310 * t283 - t307 * t284 + t306 * t287 + t308 * t288, t309 * t331 + t311 * (-t277 * qJD(3) - t283 * t304) - t323 * t273, t273, -t284, (-t273 * t300 - t284 * t303) * r_i_i_C(1) + (-t273 * t303 + t284 * t300) * r_i_i_C(2) + ((-t277 * t303 - t287 * t300) * r_i_i_C(1) + (t277 * t300 - t287 * t303) * r_i_i_C(2)) * qJD(6); 0, ((-t307 * qJD(2) + t308) * t302 + (-t310 * qJD(2) + t306) * t305) * t297, t309 * t312 + t311 * (-t291 * qJD(3) + t304 * t320) - t323 * t281, t281, -t319, (-t281 * t300 - t303 * t319) * r_i_i_C(1) + (-t281 * t303 + t300 * t319) * r_i_i_C(2) + ((-t291 * t303 - t300 * t327) * r_i_i_C(1) + (t291 * t300 - t303 * t327) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end