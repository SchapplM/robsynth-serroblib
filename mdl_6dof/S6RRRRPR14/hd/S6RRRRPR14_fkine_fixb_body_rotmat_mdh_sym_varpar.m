% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRPR14 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:41
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRRPR14_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR14_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:41:50
	% EndTime: 2020-11-04 22:41:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:41:50
	% EndTime: 2020-11-04 22:41:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t189 = cos(qJ(1));
	t188 = sin(qJ(1));
	t1 = [t189, -t188, 0, 0; t188, t189, 0, 0; 0, 0, 1, pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:41:50
	% EndTime: 2020-11-04 22:41:50
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t190 = sin(pkin(6));
	t193 = sin(qJ(1));
	t201 = t193 * t190;
	t192 = sin(qJ(2));
	t200 = t193 * t192;
	t194 = cos(qJ(2));
	t199 = t193 * t194;
	t195 = cos(qJ(1));
	t198 = t195 * t190;
	t197 = t195 * t192;
	t196 = t195 * t194;
	t191 = cos(pkin(6));
	t1 = [-t191 * t200 + t196, -t191 * t199 - t197, t201, t195 * pkin(1) + pkin(9) * t201 + 0; t191 * t197 + t199, t191 * t196 - t200, -t198, t193 * pkin(1) - pkin(9) * t198 + 0; t190 * t192, t190 * t194, t191, t191 * pkin(9) + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:41:50
	% EndTime: 2020-11-04 22:41:51
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (44->32), mult. (104->58), div. (0->0), fcn. (138->10), ass. (0->29)
	t208 = sin(pkin(7));
	t211 = cos(pkin(6));
	t230 = t208 * t211;
	t216 = cos(qJ(2));
	t229 = t208 * t216;
	t209 = sin(pkin(6));
	t210 = cos(pkin(7));
	t228 = t209 * t210;
	t227 = t211 * t216;
	t212 = sin(qJ(3));
	t213 = sin(qJ(2));
	t226 = t212 * t213;
	t225 = t212 * t216;
	t215 = cos(qJ(3));
	t224 = t213 * t215;
	t214 = sin(qJ(1));
	t223 = t214 * t213;
	t222 = t215 * t216;
	t217 = cos(qJ(1));
	t221 = t217 * t213;
	t220 = t217 * t216;
	t219 = -pkin(2) * t213 + pkin(10) * t229;
	t203 = -t208 * t209 + t210 * t227;
	t218 = t203 * t212 + t211 * t224;
	t207 = t210 * pkin(10) + pkin(9);
	t205 = t208 * t213 * pkin(10) + pkin(2) * t216 + pkin(1);
	t204 = -t210 * t226 + t222;
	t202 = t209 * t207 + t219 * t211;
	t1 = [t217 * t204 - t218 * t214, (-t203 * t214 - t210 * t221) * t215 - (-t211 * t223 + t220) * t212, (t214 * t227 + t221) * t208 + t214 * t228, t202 * t214 + t205 * t217 + 0; t214 * t204 + t218 * t217, (t203 * t215 - t211 * t226) * t217 - t214 * (t210 * t224 + t225), -(t211 * t220 - t223) * t208 - t217 * t228, -t202 * t217 + t205 * t214 + 0; t212 * t230 + (t210 * t225 + t224) * t209, t215 * t230 + (t210 * t222 - t226) * t209, -t209 * t229 + t211 * t210, t207 * t211 - t219 * t209 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:41:51
	% EndTime: 2020-11-04 22:41:51
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (91->44), mult. (209->80), div. (0->0), fcn. (264->12), ass. (0->38)
	t246 = sin(pkin(7));
	t248 = cos(pkin(7));
	t251 = sin(qJ(3));
	t255 = cos(qJ(3));
	t261 = pkin(3) * t251 - pkin(11) * t255;
	t236 = t246 * pkin(10) - t261 * t248;
	t241 = t255 * pkin(3) + pkin(11) * t251 + pkin(2);
	t252 = sin(qJ(2));
	t256 = cos(qJ(2));
	t273 = -t236 * t256 + t241 * t252;
	t247 = sin(pkin(6));
	t270 = t246 * t247;
	t269 = t246 * t251;
	t268 = t246 * t252;
	t249 = cos(pkin(6));
	t267 = t249 * t256;
	t266 = t251 * t252;
	t265 = t251 * t256;
	t264 = t252 * t255;
	t257 = cos(qJ(1));
	t263 = t252 * t257;
	t262 = t255 * t256;
	t259 = t248 * t265 + t264;
	t233 = -t247 * t269 + t259 * t249;
	t237 = t246 * t267 + t247 * t248;
	t250 = sin(qJ(4));
	t254 = cos(qJ(4));
	t260 = t233 * t250 + t254 * t237;
	t258 = t248 * pkin(10) + t261 * t246 + pkin(9);
	t253 = sin(qJ(1));
	t240 = t248 * t266 - t262;
	t239 = t248 * t267 - t270;
	t238 = -t249 * t248 + t256 * t270;
	t235 = t240 * t250 + t254 * t268;
	t234 = t259 * t247 + t249 * t269;
	t232 = t236 * t252 + t241 * t256 + pkin(1);
	t231 = t247 * t258 - t273 * t249;
	t1 = [(-t233 * t253 - t257 * t240) * t254 + (t237 * t253 + t246 * t263) * t250, t257 * t235 + t260 * t253, (t239 * t253 + t248 * t263) * t255 + (-t253 * t249 * t252 + t257 * t256) * t251, t231 * t253 + t232 * t257 + 0; (t233 * t254 - t250 * t237) * t257 + (-t240 * t254 + t250 * t268) * t253, t235 * t253 - t260 * t257, (-t239 * t255 + t249 * t266) * t257 + t253 * (t248 * t264 + t265), -t231 * t257 + t232 * t253 + 0; t234 * t254 - t250 * t238, -t234 * t250 - t254 * t238, -t249 * t246 * t255 + (-t248 * t262 + t266) * t247, t273 * t247 + t258 * t249 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:41:51
	% EndTime: 2020-11-04 22:41:51
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (167->60), mult. (379->109), div. (0->0), fcn. (468->14), ass. (0->51)
	t294 = sin(pkin(7));
	t300 = sin(qJ(3));
	t304 = cos(qJ(3));
	t299 = sin(qJ(4));
	t303 = cos(qJ(4));
	t340 = pkin(4) * t303 + qJ(5) * t299 + pkin(3);
	t342 = -pkin(11) * t304 + t300 * t340;
	t347 = -t342 * t294 - pkin(9);
	t305 = cos(qJ(2));
	t297 = cos(pkin(7));
	t324 = t300 * t303;
	t287 = -t294 * t299 + t297 * t324;
	t293 = sin(pkin(13));
	t296 = cos(pkin(13));
	t326 = t297 * t304;
	t313 = -t296 * t287 + t293 * t326;
	t346 = t313 * t305;
	t316 = t299 * pkin(4) - qJ(5) * t303 + pkin(10);
	t278 = -t294 * t316 + t297 * t342;
	t285 = pkin(11) * t300 + t304 * t340 + pkin(2);
	t301 = sin(qJ(2));
	t345 = t278 * t305 + t285 * t301;
	t295 = sin(pkin(6));
	t298 = cos(pkin(6));
	t323 = t301 * t303;
	t321 = t298 * t323;
	t325 = t300 * t301;
	t331 = t294 * t296;
	t281 = t293 * t287 + t296 * t326;
	t334 = t281 * t305;
	t327 = t297 * t299;
	t310 = t294 * t324 + t327;
	t341 = t310 * t295;
	t344 = t293 * t341 + (t296 * t325 - t334) * t298 - (t293 * t321 - t295 * t331) * t304;
	t343 = t296 * t341 - (t293 * t325 - t346) * t298 - (t295 * t294 * t293 + t296 * t321) * t304;
	t330 = t294 * t300;
	t329 = t294 * t305;
	t328 = t295 * t297;
	t322 = t303 * t304;
	t309 = t297 * t300 * t305 + t301 * t304;
	t317 = (-t295 * t330 + t309 * t298) * t299 + t303 * (t298 * t329 + t328);
	t312 = t293 * t300 + t296 * t322;
	t311 = -t293 * t322 + t296 * t300;
	t306 = cos(qJ(1));
	t302 = sin(qJ(1));
	t282 = (t297 * t325 - t304 * t305) * t299 + t294 * t323;
	t277 = -t301 * t313 - t312 * t305;
	t276 = t281 * t301 + t311 * t305;
	t275 = -t278 * t301 + t285 * t305 + pkin(1);
	t274 = t295 * t347 + t345 * t298 - t316 * t328;
	t1 = [-t277 * t306 + t343 * t302, t276 * t306 - t344 * t302, -t306 * t282 - t317 * t302, -t274 * t302 + t275 * t306 + 0; -t277 * t302 - t343 * t306, t302 * t276 + t344 * t306, -t282 * t302 + t317 * t306, t274 * t306 + t275 * t302 + 0; (t296 * t327 + (-t293 * t304 + t296 * t324) * t294) * t298 + (t312 * t301 - t346) * t295, (-t310 * t293 - t304 * t331) * t298 + (t311 * t301 - t334) * t295, (t309 * t295 + t298 * t330) * t299 + t303 * (t295 * t329 - t298 * t297), pkin(8) + 0 + t345 * t295 + (t316 * t297 - t347) * t298; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:41:51
	% EndTime: 2020-11-04 22:41:51
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (234->60), mult. (401->99), div. (0->0), fcn. (490->16), ass. (0->54)
	t367 = cos(pkin(13)) * pkin(5) + pkin(4);
	t375 = pkin(12) + qJ(5);
	t376 = sin(qJ(4));
	t380 = cos(qJ(4));
	t362 = t367 * t380 + t375 * t376 + pkin(3);
	t366 = sin(pkin(13)) * pkin(5) + pkin(11);
	t377 = sin(qJ(3));
	t381 = cos(qJ(3));
	t355 = t362 * t381 + t366 * t377 + pkin(2);
	t378 = sin(qJ(2));
	t403 = t355 * t378;
	t371 = sin(pkin(7));
	t402 = t371 * t377;
	t401 = t371 * t381;
	t382 = cos(qJ(2));
	t400 = t371 * t382;
	t399 = t377 * t380;
	t398 = t377 * t382;
	t397 = t378 * t377;
	t396 = t378 * t381;
	t395 = t381 * t382;
	t372 = sin(pkin(6));
	t374 = cos(pkin(6));
	t373 = cos(pkin(7));
	t388 = t371 * t399 + t373 * t376;
	t364 = -t371 * t376 + t373 * t399;
	t389 = t364 * t382 + t380 * t396;
	t352 = -t372 * t388 + t389 * t374;
	t359 = t378 * t364 - t380 * t395;
	t379 = sin(qJ(1));
	t383 = cos(qJ(1));
	t394 = t352 * t383 - t359 * t379;
	t393 = t352 * t379 + t383 * t359;
	t386 = t373 * t398 + t396;
	t392 = (-t372 * t402 + t386 * t374) * t376 + t380 * (t372 * t373 + t374 * t400);
	t391 = t362 * t377 - t366 * t381;
	t390 = t367 * t376 - t375 * t380 + pkin(10);
	t387 = t373 * t395 - t397;
	t385 = t391 * t373;
	t384 = t390 * t372;
	t370 = pkin(13) + qJ(6);
	t369 = cos(t370);
	t368 = sin(t370);
	t365 = t373 * t396 + t398;
	t360 = (t373 * t397 - t395) * t376 + t371 * t380 * t378;
	t358 = t387 * t372 + t374 * t401;
	t357 = -t372 * t401 + t387 * t374;
	t354 = t357 * t383 - t379 * t365;
	t353 = t357 * t379 + t383 * t365;
	t351 = t389 * t372 + t374 * t388;
	t350 = -t371 * t390 + t385;
	t349 = -t350 * t378 + t355 * t382 + pkin(1);
	t348 = (-t391 * t371 - pkin(9)) * t372 + (t350 * t382 + t403) * t374 - t373 * t384;
	t1 = [t368 * t353 - t393 * t369, t353 * t369 + t393 * t368, -t383 * t360 - t392 * t379, -t348 * t379 + t349 * t383 + 0; -t368 * t354 + t394 * t369, -t354 * t369 - t394 * t368, -t360 * t379 + t392 * t383, t348 * t383 + t349 * t379 + 0; t351 * t369 - t368 * t358, -t351 * t368 - t369 * t358, (t386 * t372 + t374 * t402) * t376 + t380 * (t372 * t400 - t374 * t373), 0 + pkin(8) + (t390 * t373 + pkin(9)) * t374 + (t382 * t385 + t403) * t372 + (t391 * t374 - t382 * t384) * t371; 0, 0, 0, 1;];
	Tc_mdh = t1;
end