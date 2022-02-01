% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPRR6
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
%   Siehe auch: S5RRPRR6_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRPRR6_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR6_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
JRD_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:18:05
	% EndTime: 2022-01-20 11:18:05
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:18:05
	% EndTime: 2022-01-20 11:18:06
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:18:06
	% EndTime: 2022-01-20 11:18:06
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->8), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t47 = qJD(1) + qJD(2);
	t48 = qJ(1) + qJ(2);
	t49 = t47 * cos(t48);
	t44 = t47 * sin(t48);
	t1 = [-t49, -t49, 0, 0, 0; -t44, -t44, 0, 0, 0; 0, 0, 0, 0, 0; t44, t44, 0, 0, 0; -t49, -t49, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:18:05
	% EndTime: 2022-01-20 11:18:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (30->8), mult. (20->8), div. (0->0), fcn. (20->4), ass. (0->13)
	t55 = qJ(1) + qJ(2);
	t52 = sin(t55);
	t54 = qJD(1) + qJD(2);
	t62 = t54 * t52;
	t61 = t54 * sin(pkin(9));
	t60 = t54 * cos(pkin(9));
	t59 = t52 * t60;
	t53 = cos(t55);
	t58 = t53 * t60;
	t51 = t54 * t53;
	t50 = t53 * t61;
	t49 = t52 * t61;
	t1 = [-t58, -t58, 0, 0, 0; -t59, -t59, 0, 0, 0; 0, 0, 0, 0, 0; t50, t50, 0, 0, 0; t49, t49, 0, 0, 0; 0, 0, 0, 0, 0; -t62, -t62, 0, 0, 0; t51, t51, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:18:06
	% EndTime: 2022-01-20 11:18:06
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (121->20), mult. (132->28), div. (0->0), fcn. (132->6), ass. (0->24)
	t199 = qJ(1) + qJ(2);
	t196 = sin(t199);
	t203 = cos(qJ(4));
	t214 = t196 * t203;
	t197 = cos(t199);
	t202 = sin(qJ(4));
	t213 = t197 * t202;
	t198 = qJD(1) + qJD(2);
	t200 = sin(pkin(9));
	t212 = t198 * t200;
	t201 = cos(pkin(9));
	t211 = t201 * t202;
	t210 = t201 * t203;
	t209 = qJD(4) * t202;
	t208 = qJD(4) * t203;
	t207 = t196 * t212;
	t206 = t197 * t212;
	t205 = -t196 * t202 - t197 * t210;
	t204 = t196 * t211 + t197 * t203;
	t193 = t204 * qJD(4) + t205 * t198;
	t192 = (t197 * t211 - t214) * t198 + (t196 * t210 - t213) * qJD(4);
	t191 = -t198 * t213 - t196 * t208 + (t197 * t209 + t198 * t214) * t201;
	t190 = t205 * qJD(4) + t204 * t198;
	t1 = [t193, t193, 0, t190, 0; -t191, -t191, 0, -t192, 0; 0, 0, 0, -t200 * t208, 0; t192, t192, 0, t191, 0; t190, t190, 0, t193, 0; 0, 0, 0, t200 * t209, 0; -t206, -t206, 0, 0, 0; -t207, -t207, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:18:06
	% EndTime: 2022-01-20 11:18:06
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (262->22), mult. (176->22), div. (0->0), fcn. (176->6), ass. (0->29)
	t267 = qJ(4) + qJ(5);
	t261 = sin(t267);
	t268 = qJ(1) + qJ(2);
	t262 = sin(t268);
	t283 = t261 * t262;
	t264 = cos(t268);
	t282 = t261 * t264;
	t263 = cos(t267);
	t281 = t262 * t263;
	t280 = t263 * t264;
	t265 = qJD(4) + qJD(5);
	t269 = sin(pkin(9));
	t279 = t265 * t269;
	t266 = qJD(1) + qJD(2);
	t278 = t266 * t269;
	t277 = t262 * t278;
	t276 = t264 * t278;
	t275 = t263 * t279;
	t270 = cos(pkin(9));
	t274 = t266 * t270 - t265;
	t273 = t265 * t270 - t266;
	t272 = t265 * t281 + t266 * t282;
	t271 = t265 * t282 + t266 * t281;
	t260 = t261 * t279;
	t255 = t273 * t283 - t274 * t280;
	t254 = t272 * t270 - t271;
	t253 = t271 * t270 - t272;
	t252 = -t273 * t280 + t274 * t283;
	t1 = [t255, t255, 0, t252, t252; -t253, -t253, 0, -t254, -t254; 0, 0, 0, -t275, -t275; t254, t254, 0, t253, t253; t252, t252, 0, t255, t255; 0, 0, 0, t260, t260; -t276, -t276, 0, 0, 0; -t277, -t277, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
end