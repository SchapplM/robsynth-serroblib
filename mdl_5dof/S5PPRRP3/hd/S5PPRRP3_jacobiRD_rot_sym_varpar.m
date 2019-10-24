% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PPRRP3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:20
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PPRRP3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP3_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRP3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:04
	% EndTime: 2019-10-24 10:20:04
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:04
	% EndTime: 2019-10-24 10:20:04
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:04
	% EndTime: 2019-10-24 10:20:04
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:04
	% EndTime: 2019-10-24 10:20:04
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (24->15), div. (0->0), fcn. (24->6), ass. (0->11)
	t66 = sin(pkin(7));
	t69 = sin(qJ(3));
	t75 = t66 * t69;
	t70 = cos(qJ(3));
	t74 = t66 * t70;
	t68 = cos(pkin(7));
	t73 = t68 * t69;
	t72 = t68 * t70;
	t71 = qJD(3) * sin(pkin(8));
	t67 = cos(pkin(8));
	t1 = [0, 0, (-t67 * t72 - t75) * qJD(3), 0, 0; 0, 0, (-t67 * t74 + t73) * qJD(3), 0, 0; 0, 0, -t70 * t71, 0, 0; 0, 0, (t67 * t73 - t74) * qJD(3), 0, 0; 0, 0, (t67 * t75 + t72) * qJD(3), 0, 0; 0, 0, t69 * t71, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:04
	% EndTime: 2019-10-24 10:20:04
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (37->23), mult. (140->57), div. (0->0), fcn. (148->8), ass. (0->28)
	t218 = sin(pkin(8));
	t222 = sin(qJ(4));
	t236 = t218 * t222;
	t224 = cos(qJ(4));
	t235 = t218 * t224;
	t225 = cos(qJ(3));
	t234 = t218 * t225;
	t219 = sin(pkin(7));
	t223 = sin(qJ(3));
	t233 = t219 * t223;
	t232 = t219 * t225;
	t221 = cos(pkin(7));
	t231 = t221 * t223;
	t230 = t221 * t225;
	t229 = qJD(3) * t225;
	t228 = qJD(4) * t222;
	t227 = qJD(4) * t224;
	t226 = t218 * qJD(3) * t223;
	t220 = cos(pkin(8));
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
	% StartTime: 2019-10-24 10:20:05
	% EndTime: 2019-10-24 10:20:05
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (37->23), mult. (140->57), div. (0->0), fcn. (148->8), ass. (0->28)
	t277 = sin(pkin(8));
	t281 = sin(qJ(4));
	t295 = t277 * t281;
	t283 = cos(qJ(4));
	t294 = t277 * t283;
	t284 = cos(qJ(3));
	t293 = t277 * t284;
	t278 = sin(pkin(7));
	t282 = sin(qJ(3));
	t292 = t278 * t282;
	t291 = t278 * t284;
	t280 = cos(pkin(7));
	t290 = t280 * t282;
	t289 = t280 * t284;
	t288 = qJD(3) * t284;
	t287 = qJD(4) * t281;
	t286 = qJD(4) * t283;
	t285 = t277 * qJD(3) * t282;
	t279 = cos(pkin(8));
	t276 = t279 * t289 + t292;
	t275 = -t279 * t290 + t291;
	t274 = t279 * t291 - t290;
	t273 = -t279 * t292 - t289;
	t272 = t276 * qJD(3);
	t271 = t275 * qJD(3);
	t270 = t274 * qJD(3);
	t269 = t273 * qJD(3);
	t1 = [0, 0, -t272 * t283 - t275 * t287, -t271 * t281 + (-t276 * t283 - t280 * t295) * qJD(4), 0; 0, 0, -t270 * t283 - t273 * t287, -t269 * t281 + (-t274 * t283 - t278 * t295) * qJD(4), 0; 0, 0, (t282 * t287 - t283 * t288) * t277, t281 * t285 + (t279 * t281 - t283 * t293) * qJD(4), 0; 0, 0, t271, 0, 0; 0, 0, t269, 0, 0; 0, 0, -t285, 0, 0; 0, 0, -t272 * t281 + t275 * t286, t271 * t283 + (-t276 * t281 + t280 * t294) * qJD(4), 0; 0, 0, -t270 * t281 + t273 * t286, t269 * t283 + (-t274 * t281 + t278 * t294) * qJD(4), 0; 0, 0, (-t281 * t288 - t282 * t286) * t277, -t283 * t285 + (-t279 * t283 - t281 * t293) * qJD(4), 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end