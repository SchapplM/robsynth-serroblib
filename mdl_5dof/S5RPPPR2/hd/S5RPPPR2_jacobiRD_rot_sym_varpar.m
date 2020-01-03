% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPPPR2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:23:58
	% EndTime: 2020-01-03 11:23:58
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:23:58
	% EndTime: 2020-01-03 11:23:58
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t10 = qJD(1) * sin(qJ(1));
	t9 = qJD(1) * cos(qJ(1));
	t1 = [0, 0, 0, 0, 0; -t10, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0; -t10, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:23:58
	% EndTime: 2020-01-03 11:23:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t43 = qJD(1) * sin(qJ(1));
	t42 = qJD(1) * cos(qJ(1));
	t39 = cos(pkin(7));
	t38 = sin(pkin(7));
	t1 = [0, 0, 0, 0, 0; -t39 * t43, 0, 0, 0, 0; t39 * t42, 0, 0, 0, 0; 0, 0, 0, 0, 0; t38 * t43, 0, 0, 0, 0; -t38 * t42, 0, 0, 0, 0; 0, 0, 0, 0, 0; t42, 0, 0, 0, 0; t43, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:23:58
	% EndTime: 2020-01-03 11:23:58
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t87 = cos(pkin(7));
	t88 = sin(qJ(1));
	t92 = t87 * t88;
	t89 = cos(qJ(1));
	t91 = t87 * t89;
	t90 = qJD(1) * sin(pkin(7));
	t86 = cos(pkin(8));
	t84 = sin(pkin(8));
	t1 = [0, 0, 0, 0, 0; (t84 * t89 - t86 * t92) * qJD(1), 0, 0, 0, 0; (t84 * t88 + t86 * t91) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0; (t84 * t92 + t86 * t89) * qJD(1), 0, 0, 0, 0; (-t84 * t91 + t86 * t88) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0; -t88 * t90, 0, 0, 0, 0; t89 * t90, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:23:58
	% EndTime: 2020-01-03 11:23:58
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (10->8), mult. (46->25), div. (0->0), fcn. (46->8), ass. (0->15)
	t144 = cos(pkin(7));
	t145 = sin(qJ(1));
	t151 = t144 * t145;
	t146 = cos(qJ(1));
	t150 = t144 * t146;
	t149 = qJD(1) * sin(pkin(7));
	t148 = t145 * t149;
	t147 = t146 * t149;
	t143 = cos(pkin(8));
	t142 = cos(pkin(9));
	t140 = sin(pkin(8));
	t139 = sin(pkin(9));
	t138 = (t140 * t145 + t143 * t150) * qJD(1);
	t137 = (t140 * t146 - t143 * t151) * qJD(1);
	t1 = [0, 0, 0, 0, 0; t137 * t142 - t139 * t148, 0, 0, 0, 0; t138 * t142 + t139 * t147, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t137 * t139 - t142 * t148, 0, 0, 0, 0; -t138 * t139 + t142 * t147, 0, 0, 0, 0; 0, 0, 0, 0, 0; (-t140 * t151 - t143 * t146) * qJD(1), 0, 0, 0, 0; (t140 * t150 - t143 * t145) * qJD(1), 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:23:59
	% EndTime: 2020-01-03 11:23:59
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (80->28), mult. (266->60), div. (0->0), fcn. (292->10), ass. (0->37)
	t281 = sin(pkin(8));
	t285 = cos(pkin(7));
	t284 = cos(pkin(8));
	t287 = sin(qJ(1));
	t296 = t287 * t284;
	t292 = t285 * t296;
	t289 = cos(qJ(1));
	t293 = qJD(1) * t289;
	t271 = qJD(1) * t292 - t281 * t293;
	t280 = sin(pkin(9));
	t283 = cos(pkin(9));
	t282 = sin(pkin(7));
	t298 = t282 * t287;
	t291 = qJD(1) * t298;
	t264 = t271 * t283 + t280 * t291;
	t294 = t289 * t284;
	t297 = t287 * t281;
	t277 = t285 * t294 + t297;
	t268 = t289 * t282 * t280 + t277 * t283;
	t275 = t285 * t297 + t294;
	t270 = t275 * qJD(1);
	t295 = t289 * t281;
	t276 = t285 * t295 - t296;
	t286 = sin(qJ(5));
	t288 = cos(qJ(5));
	t307 = (t268 * t288 + t276 * t286) * qJD(5) - t264 * t286 + t270 * t288;
	t304 = -t264 * t288 - t270 * t286 + (-t268 * t286 + t276 * t288) * qJD(5);
	t299 = t281 * t282;
	t290 = t282 * t293;
	t274 = t282 * t284 * t283 - t285 * t280;
	t273 = t277 * qJD(1);
	t272 = t276 * qJD(1);
	t267 = (t292 - t295) * t283 + t280 * t298;
	t266 = t273 * t283 + t280 * t290;
	t263 = t266 * t288 + t272 * t286 + (-t267 * t286 + t275 * t288) * qJD(5);
	t262 = -t266 * t286 + t272 * t288 + (-t267 * t288 - t275 * t286) * qJD(5);
	t1 = [0, 0, 0, 0, (-t274 * t288 - t286 * t299) * qJD(5); t304, 0, 0, 0, t262; t263, 0, 0, 0, t307; 0, 0, 0, 0, (t274 * t286 - t288 * t299) * qJD(5); -t307, 0, 0, 0, -t263; t262, 0, 0, 0, t304; 0, 0, 0, 0, 0; -t271 * t280 + t283 * t291, 0, 0, 0, 0; t273 * t280 - t283 * t290, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end