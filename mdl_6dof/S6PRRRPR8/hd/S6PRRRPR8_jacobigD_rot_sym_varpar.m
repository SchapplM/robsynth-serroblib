% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRPR8
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:00
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRPR8_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR8_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR8_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_jacobigD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:04
	% EndTime: 2019-10-09 23:00:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:04
	% EndTime: 2019-10-09 23:00:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:04
	% EndTime: 2019-10-09 23:00:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:05
	% EndTime: 2019-10-09 23:00:05
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->4), mult. (15->10), div. (0->0), fcn. (15->7), ass. (0->7)
	t106 = sin(qJ(2));
	t109 = cos(pkin(6)) * t106;
	t108 = qJD(2) * sin(pkin(7));
	t107 = cos(qJ(2));
	t104 = cos(pkin(12));
	t102 = sin(pkin(12));
	t1 = [0, 0, -(t102 * t109 - t104 * t107) * t108, 0, 0, 0; 0, 0, -(-t102 * t107 - t104 * t109) * t108, 0, 0, 0; 0, 0, sin(pkin(6)) * t106 * t108, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:05
	% EndTime: 2019-10-09 23:00:05
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (24->17), mult. (88->43), div. (0->0), fcn. (92->10), ass. (0->24)
	t182 = sin(pkin(7));
	t183 = sin(pkin(6));
	t201 = t183 * t182;
	t185 = cos(pkin(7));
	t189 = cos(qJ(3));
	t200 = t185 * t189;
	t186 = cos(pkin(6));
	t188 = sin(qJ(2));
	t199 = t186 * t188;
	t190 = cos(qJ(2));
	t198 = t186 * t190;
	t187 = sin(qJ(3));
	t197 = t187 * t190;
	t196 = t188 * t189;
	t195 = qJD(2) * t187;
	t181 = sin(pkin(12));
	t184 = cos(pkin(12));
	t194 = -t181 * t188 + t184 * t198;
	t193 = t181 * t190 + t184 * t199;
	t192 = -t181 * t198 - t184 * t188;
	t191 = t181 * t199 - t184 * t190;
	t180 = t191 * qJD(2);
	t179 = t193 * qJD(2);
	t1 = [0, 0, -t180 * t182, -t180 * t200 + t192 * t195 + (-t191 * t189 + (t181 * t201 + t185 * t192) * t187) * qJD(3), 0, 0; 0, 0, t179 * t182, t179 * t200 + t194 * t195 + (t193 * t189 + (-t184 * t201 + t185 * t194) * t187) * qJD(3), 0, 0; 0, 0, qJD(2) * t188 * t201, t186 * t182 * qJD(3) * t187 + ((t185 * t197 + t196) * qJD(3) + (t185 * t196 + t197) * qJD(2)) * t183, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:05
	% EndTime: 2019-10-09 23:00:05
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (24->17), mult. (88->43), div. (0->0), fcn. (92->10), ass. (0->24)
	t203 = sin(pkin(7));
	t204 = sin(pkin(6));
	t222 = t204 * t203;
	t206 = cos(pkin(7));
	t210 = cos(qJ(3));
	t221 = t206 * t210;
	t207 = cos(pkin(6));
	t209 = sin(qJ(2));
	t220 = t207 * t209;
	t211 = cos(qJ(2));
	t219 = t207 * t211;
	t208 = sin(qJ(3));
	t218 = t208 * t211;
	t217 = t209 * t210;
	t216 = qJD(2) * t208;
	t202 = sin(pkin(12));
	t205 = cos(pkin(12));
	t215 = -t202 * t209 + t205 * t219;
	t214 = t202 * t211 + t205 * t220;
	t213 = -t202 * t219 - t205 * t209;
	t212 = t202 * t220 - t205 * t211;
	t201 = t212 * qJD(2);
	t200 = t214 * qJD(2);
	t1 = [0, 0, -t201 * t203, -t201 * t221 + t213 * t216 + (-t212 * t210 + (t202 * t222 + t213 * t206) * t208) * qJD(3), 0, 0; 0, 0, t200 * t203, t200 * t221 + t215 * t216 + (t214 * t210 + (-t205 * t222 + t215 * t206) * t208) * qJD(3), 0, 0; 0, 0, qJD(2) * t209 * t222, t207 * t203 * qJD(3) * t208 + ((t206 * t218 + t217) * qJD(3) + (t206 * t217 + t218) * qJD(2)) * t204, 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:06
	% EndTime: 2019-10-09 23:00:06
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (70->40), mult. (240->92), div. (0->0), fcn. (263->12), ass. (0->40)
	t259 = sin(pkin(7));
	t260 = sin(pkin(6));
	t289 = t259 * t260;
	t264 = sin(qJ(4));
	t288 = t259 * t264;
	t265 = sin(qJ(3));
	t287 = t259 * t265;
	t266 = sin(qJ(2));
	t286 = t259 * t266;
	t262 = cos(pkin(7));
	t285 = t260 * t262;
	t284 = t262 * t265;
	t268 = cos(qJ(3));
	t283 = t262 * t268;
	t263 = cos(pkin(6));
	t282 = t263 * t266;
	t269 = cos(qJ(2));
	t281 = t263 * t269;
	t280 = t265 * t266;
	t279 = t265 * t269;
	t278 = t266 * t268;
	t277 = t268 * t269;
	t267 = cos(qJ(4));
	t276 = qJD(3) * t267;
	t258 = sin(pkin(12));
	t261 = cos(pkin(12));
	t254 = -t258 * t266 + t261 * t281;
	t275 = t254 * t262 - t261 * t289;
	t256 = -t258 * t281 - t261 * t266;
	t274 = t256 * t262 + t258 * t289;
	t255 = t258 * t269 + t261 * t282;
	t273 = t258 * t282 - t261 * t269;
	t272 = t262 * t279 + t278;
	t271 = t255 * t268 + t275 * t265;
	t270 = t274 * t265 - t268 * t273;
	t253 = t273 * qJD(2);
	t252 = t256 * qJD(2);
	t251 = t255 * qJD(2);
	t250 = t254 * qJD(2);
	t1 = [0, 0, -t253 * t259, t270 * qJD(3) + t252 * t265 - t253 * t283, 0, (t252 * t268 + t253 * t284) * t267 - t253 * t288 + (-t270 * t264 + (-t256 * t259 + t258 * t285) * t267) * qJD(4) + (t265 * t273 + t274 * t268) * t276; 0, 0, t251 * t259, t271 * qJD(3) + t250 * t265 + t251 * t283, 0, (t250 * t268 - t251 * t284) * t267 + t251 * t288 + (-t271 * t264 + (-t254 * t259 - t261 * t285) * t267) * qJD(4) + (-t255 * t265 + t275 * t268) * t276; 0, 0, t260 * qJD(2) * t286, t263 * qJD(3) * t287 + (t272 * qJD(3) + (t262 * t278 + t279) * qJD(2)) * t260, 0, (t259 * t268 * t276 + (t262 * t267 - t264 * t287) * qJD(4)) * t263 + ((-t259 * t269 * t267 - t272 * t264) * qJD(4) + (t262 * t277 - t280) * t276 + ((-t262 * t280 + t277) * t267 + t264 * t286) * qJD(2)) * t260;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end