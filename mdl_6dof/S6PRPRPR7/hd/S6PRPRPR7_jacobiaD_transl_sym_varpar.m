% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:40
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRPR7_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR7_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR7_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.05s
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
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (24->13), mult. (82->26), div. (0->0), fcn. (72->6), ass. (0->16)
	t135 = -pkin(2) + r_i_i_C(2);
	t134 = r_i_i_C(3) + qJ(3);
	t126 = cos(pkin(6));
	t127 = sin(qJ(2));
	t133 = t126 * t127;
	t128 = cos(qJ(2));
	t132 = t126 * t128;
	t131 = qJD(2) * t134;
	t123 = sin(pkin(10));
	t125 = cos(pkin(10));
	t130 = t123 * t128 + t125 * t133;
	t129 = t123 * t133 - t125 * t128;
	t124 = sin(pkin(6));
	t122 = t129 * qJD(2);
	t120 = t130 * qJD(2);
	t1 = [0, -t129 * qJD(3) - t135 * t122 - (t123 * t132 + t125 * t127) * t131, -t122, 0, 0, 0; 0, t130 * qJD(3) + t135 * t120 - (t123 * t127 - t125 * t132) * t131, t120, 0, 0, 0; 0, (qJD(3) * t127 + (t135 * t127 + t134 * t128) * qJD(2)) * t124, t124 * qJD(2) * t127, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:57
	% EndTime: 2019-10-09 21:40:57
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (67->29), mult. (228->61), div. (0->0), fcn. (208->8), ass. (0->28)
	t178 = sin(pkin(6));
	t181 = sin(qJ(4));
	t199 = t178 * t181;
	t183 = cos(qJ(4));
	t198 = t178 * t183;
	t184 = cos(qJ(2));
	t197 = t178 * t184;
	t179 = cos(pkin(10));
	t196 = t179 * t184;
	t180 = cos(pkin(6));
	t182 = sin(qJ(2));
	t195 = t180 * t182;
	t194 = t180 * t184;
	t193 = qJD(2) * t182;
	t192 = -pkin(2) - pkin(8) - r_i_i_C(3);
	t177 = sin(pkin(10));
	t191 = t177 * t193;
	t190 = t178 * t193;
	t189 = qJD(2) * t196;
	t188 = t183 * r_i_i_C(1) - t181 * r_i_i_C(2);
	t187 = t181 * r_i_i_C(1) + t183 * r_i_i_C(2) + qJ(3);
	t186 = t177 * t184 + t179 * t195;
	t173 = t177 * t194 + t179 * t182;
	t185 = t188 * qJD(4) + qJD(3);
	t171 = t177 * t182 - t179 * t194;
	t170 = -t180 * t191 + t189;
	t168 = t186 * qJD(2);
	t1 = [0, t185 * (-t177 * t195 + t196) + t192 * t170 - t187 * t173 * qJD(2), t170, t188 * t170 + ((-t173 * t181 - t177 * t198) * r_i_i_C(1) + (-t173 * t183 + t177 * t199) * r_i_i_C(2)) * qJD(4), 0, 0; 0, t185 * t186 + t192 * t168 - t187 * (-t180 * t189 + t191), t168, t188 * t168 + ((-t171 * t181 + t179 * t198) * r_i_i_C(1) + (-t171 * t183 - t179 * t199) * r_i_i_C(2)) * qJD(4), 0, 0; 0, (t185 * t182 + (t192 * t182 + t187 * t184) * qJD(2)) * t178, t190, t188 * t190 + ((-t180 * t183 + t181 * t197) * r_i_i_C(1) + (t180 * t181 + t183 * t197) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:57
	% EndTime: 2019-10-09 21:40:57
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (136->40), mult. (441->73), div. (0->0), fcn. (420->8), ass. (0->35)
	t226 = sin(qJ(4));
	t228 = cos(qJ(4));
	t247 = r_i_i_C(3) + qJ(5);
	t248 = -r_i_i_C(2) + pkin(4);
	t249 = t248 * t226 - t247 * t228 + qJ(3);
	t222 = sin(pkin(10));
	t224 = cos(pkin(10));
	t227 = sin(qJ(2));
	t225 = cos(pkin(6));
	t229 = cos(qJ(2));
	t241 = t225 * t229;
	t216 = t222 * t227 - t224 * t241;
	t246 = t216 * t226;
	t223 = sin(pkin(6));
	t245 = t223 * t226;
	t244 = t223 * t228;
	t243 = t223 * t229;
	t242 = t225 * t227;
	t240 = qJD(2) * t227;
	t239 = qJD(2) * t229;
	t238 = -pkin(2) - pkin(8) - r_i_i_C(1);
	t237 = t222 * t240;
	t236 = t223 * t240;
	t235 = t224 * t239;
	t218 = t222 * t241 + t224 * t227;
	t234 = t218 * t226 + t222 * t244;
	t233 = t222 * t229 + t224 * t242;
	t232 = -t225 * t228 + t226 * t243;
	t230 = -qJD(5) * t228 + qJD(3) + (t247 * t226 + t248 * t228) * qJD(4);
	t215 = -t225 * t237 + t235;
	t213 = t233 * qJD(2);
	t210 = -t232 * qJD(4) - t228 * t236;
	t208 = -qJD(4) * t246 + (qJD(4) * t223 * t224 + t213) * t228;
	t206 = t234 * qJD(4) - t215 * t228;
	t1 = [0, t238 * t215 - t249 * t218 * qJD(2) + t230 * (-t222 * t242 + t224 * t229), t215, t234 * qJD(5) - t248 * t206 - t247 * (-t215 * t226 + (-t218 * t228 + t222 * t245) * qJD(4)), t206, 0; 0, t238 * t213 - t249 * (-t225 * t235 + t237) + t230 * t233, t213, -(t224 * t244 - t246) * qJD(5) + t248 * t208 + t247 * (t213 * t226 + (t216 * t228 + t224 * t245) * qJD(4)), -t208, 0; 0, (t249 * t239 + (t238 * qJD(2) + t230) * t227) * t223, t236, -t232 * qJD(5) - t248 * t210 - t247 * (-t226 * t236 + (t225 * t226 + t228 * t243) * qJD(4)), t210, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:58
	% EndTime: 2019-10-09 21:40:58
	% DurationCPUTime: 0.49s
	% Computational Cost: add. (287->69), mult. (917->127), div. (0->0), fcn. (916->10), ass. (0->51)
	t297 = sin(qJ(6));
	t300 = cos(qJ(6));
	t312 = t297 * r_i_i_C(1) + t300 * r_i_i_C(2);
	t308 = qJD(6) * t312;
	t313 = t300 * r_i_i_C(1) - t297 * r_i_i_C(2);
	t305 = t313 * qJD(6) + qJD(5);
	t298 = sin(qJ(4));
	t301 = cos(qJ(4));
	t311 = qJ(5) + t312;
	t322 = pkin(4) + pkin(9) + r_i_i_C(3);
	t332 = t322 * t298 - t311 * t301 + qJ(3);
	t331 = cos(pkin(6));
	t294 = sin(pkin(10));
	t296 = cos(pkin(10));
	t299 = sin(qJ(2));
	t302 = cos(qJ(2));
	t317 = t302 * t331;
	t284 = t294 * t317 + t296 * t299;
	t330 = t284 * t301;
	t295 = sin(pkin(6));
	t329 = t295 * t298;
	t328 = t295 * t299;
	t327 = t295 * t301;
	t326 = t295 * t302;
	t325 = qJD(2) * t299;
	t324 = qJD(2) * t302;
	t323 = qJD(4) * t298;
	t321 = t295 * t324;
	t320 = t295 * t323;
	t319 = t295 * t325;
	t318 = t299 * t331;
	t316 = t331 * t301;
	t315 = t294 * t318;
	t314 = t296 * t317;
	t282 = t294 * t299 - t314;
	t310 = t282 * t301 + t296 * t329;
	t309 = t284 * t298 + t294 * t327;
	t307 = -pkin(2) - pkin(5) - pkin(8) - t313;
	t283 = t294 * t302 + t296 * t318;
	t286 = t331 * t298 + t301 * t326;
	t303 = qJD(3) - t305 * t301 + (t311 * t298 + t322 * t301) * qJD(4);
	t285 = t296 * t302 - t315;
	t281 = -qJD(2) * t315 + t296 * t324;
	t280 = t284 * qJD(2);
	t279 = t283 * qJD(2);
	t278 = -qJD(2) * t314 + t294 * t325;
	t272 = qJD(4) * t316 - t301 * t319 - t302 * t320;
	t267 = t294 * t329 - t330;
	t266 = -t282 * t323 + (qJD(4) * t295 * t296 + t279) * t301;
	t264 = t309 * qJD(4) - t281 * t301;
	t1 = [0, -t280 * t332 + t307 * t281 + t284 * t308 + t303 * t285, t281, t305 * t309 - t322 * t264 - t311 * (-qJD(4) * t330 - t281 * t298 + t294 * t320), t264, (t264 * t300 + t280 * t297) * r_i_i_C(1) + (-t264 * t297 + t280 * t300) * r_i_i_C(2) + ((-t267 * t297 - t285 * t300) * r_i_i_C(1) + (-t267 * t300 + t285 * t297) * r_i_i_C(2)) * qJD(6); 0, -t278 * t332 + t307 * t279 + t282 * t308 + t303 * t283, t279, -t305 * (-t282 * t298 + t296 * t327) + t322 * t266 + t311 * (t310 * qJD(4) + t279 * t298), -t266, (-t266 * t300 + t278 * t297) * r_i_i_C(1) + (t266 * t297 + t278 * t300) * r_i_i_C(2) + ((-t283 * t300 + t297 * t310) * r_i_i_C(1) + (t283 * t297 + t300 * t310) * r_i_i_C(2)) * qJD(6); 0, ((t332 * qJD(2) - t308) * t302 + (t307 * qJD(2) + t303) * t299) * t295, t319, t305 * (-t298 * t326 + t316) - t322 * t272 - t311 * (t286 * qJD(4) - t298 * t319), t272, (t272 * t300 - t297 * t321) * r_i_i_C(1) + (-t272 * t297 - t300 * t321) * r_i_i_C(2) + ((-t286 * t297 - t300 * t328) * r_i_i_C(1) + (-t286 * t300 + t297 * t328) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end