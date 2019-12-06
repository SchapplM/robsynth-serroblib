% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PPRRR3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR3_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRR3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:17:26
	% EndTime: 2019-12-05 15:17:26
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:17:26
	% EndTime: 2019-12-05 15:17:26
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:17:26
	% EndTime: 2019-12-05 15:17:27
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:17:27
	% EndTime: 2019-12-05 15:17:27
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (24->15), div. (0->0), fcn. (24->6), ass. (0->11)
	t66 = sin(pkin(8));
	t69 = sin(qJ(3));
	t75 = t66 * t69;
	t70 = cos(qJ(3));
	t74 = t66 * t70;
	t68 = cos(pkin(8));
	t73 = t68 * t69;
	t72 = t68 * t70;
	t71 = qJD(3) * sin(pkin(9));
	t67 = cos(pkin(9));
	t1 = [0, 0, (-t67 * t72 - t75) * qJD(3), 0, 0; 0, 0, (-t67 * t74 + t73) * qJD(3), 0, 0; 0, 0, -t70 * t71, 0, 0; 0, 0, (t67 * t73 - t74) * qJD(3), 0, 0; 0, 0, (t67 * t75 + t72) * qJD(3), 0, 0; 0, 0, t69 * t71, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:17:27
	% EndTime: 2019-12-05 15:17:27
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (37->23), mult. (140->57), div. (0->0), fcn. (148->8), ass. (0->28)
	t218 = sin(pkin(9));
	t222 = sin(qJ(4));
	t236 = t218 * t222;
	t224 = cos(qJ(4));
	t235 = t218 * t224;
	t225 = cos(qJ(3));
	t234 = t218 * t225;
	t219 = sin(pkin(8));
	t223 = sin(qJ(3));
	t233 = t219 * t223;
	t232 = t219 * t225;
	t221 = cos(pkin(8));
	t231 = t221 * t223;
	t230 = t221 * t225;
	t229 = qJD(3) * t225;
	t228 = qJD(4) * t222;
	t227 = qJD(4) * t224;
	t226 = t218 * qJD(3) * t223;
	t220 = cos(pkin(9));
	t217 = t220 * t230 + t233;
	t216 = -t220 * t231 + t232;
	t215 = t220 * t232 - t231;
	t214 = -t220 * t233 - t230;
	t213 = t217 * qJD(3);
	t212 = t216 * qJD(3);
	t211 = t215 * qJD(3);
	t210 = t214 * qJD(3);
	t1 = [0, 0, -t213 * t224 - t216 * t228, -t212 * t222 + (-t217 * t224 - t221 * t236) * qJD(4), 0; 0, 0, -t211 * t224 - t214 * t228, -t210 * t222 + (-t215 * t224 - t219 * t236) * qJD(4), 0; 0, 0, (t223 * t228 - t224 * t229) * t218, t222 * t226 + (t220 * t222 - t224 * t234) * qJD(4), 0; 0, 0, t213 * t222 - t216 * t227, -t212 * t224 + (t217 * t222 - t221 * t235) * qJD(4), 0; 0, 0, t211 * t222 - t214 * t227, -t210 * t224 + (t215 * t222 - t219 * t235) * qJD(4), 0; 0, 0, (t222 * t229 + t223 * t227) * t218, t224 * t226 + (t220 * t224 + t222 * t234) * qJD(4), 0; 0, 0, t212, 0, 0; 0, 0, t210, 0, 0; 0, 0, -t226, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:17:27
	% EndTime: 2019-12-05 15:17:27
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (135->22), mult. (212->48), div. (0->0), fcn. (224->8), ass. (0->39)
	t285 = qJ(4) + qJ(5);
	t282 = sin(t285);
	t284 = qJD(4) + qJD(5);
	t305 = t282 * t284;
	t283 = cos(t285);
	t304 = t283 * t284;
	t286 = sin(pkin(9));
	t303 = t284 * t286;
	t290 = sin(qJ(3));
	t302 = t284 * t290;
	t287 = sin(pkin(8));
	t301 = t287 * t290;
	t291 = cos(qJ(3));
	t300 = t287 * t291;
	t289 = cos(pkin(8));
	t299 = t289 * t290;
	t298 = t289 * t291;
	t297 = qJD(3) * t291;
	t296 = t291 * t303;
	t295 = t286 * qJD(3) * t290;
	t288 = cos(pkin(9));
	t278 = -t288 * t301 - t298;
	t274 = t278 * qJD(3);
	t294 = -t287 * t303 - t274;
	t280 = -t288 * t299 + t300;
	t276 = t280 * qJD(3);
	t293 = -t289 * t303 - t276;
	t281 = t288 * t298 + t301;
	t279 = t288 * t300 - t299;
	t292 = t284 * t288 + t295;
	t277 = t281 * qJD(3);
	t275 = t279 * qJD(3);
	t273 = t282 * t296 + t292 * t283;
	t272 = t292 * t282 - t283 * t296;
	t271 = t281 * t305 + t293 * t283;
	t270 = -t281 * t304 + t293 * t282;
	t269 = t279 * t305 + t294 * t283;
	t268 = -t279 * t304 + t294 * t282;
	t1 = [0, 0, -t277 * t283 - t280 * t305, t270, t270; 0, 0, -t275 * t283 - t278 * t305, t268, t268; 0, 0, (t282 * t302 - t283 * t297) * t286, t272, t272; 0, 0, t277 * t282 - t280 * t304, t271, t271; 0, 0, t275 * t282 - t278 * t304, t269, t269; 0, 0, (t282 * t297 + t283 * t302) * t286, t273, t273; 0, 0, t276, 0, 0; 0, 0, t274, 0, 0; 0, 0, -t295, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end