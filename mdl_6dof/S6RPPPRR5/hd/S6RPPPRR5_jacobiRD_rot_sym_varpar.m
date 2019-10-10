% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:32
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPPPRR5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPPRR5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:32:36
	% EndTime: 2019-10-09 23:32:36
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:32:36
	% EndTime: 2019-10-09 23:32:36
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:32:36
	% EndTime: 2019-10-09 23:32:36
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (1->1), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t13 = qJD(1) * sin(qJ(1));
	t11 = qJD(1) * cos(qJ(1));
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; t13, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t13, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:32:36
	% EndTime: 2019-10-09 23:32:36
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t16 = qJD(1) * sin(qJ(1));
	t15 = qJD(1) * cos(qJ(1));
	t1 = [-t16, 0, 0, 0, 0, 0; t15, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t15, 0, 0, 0, 0, 0; -t16, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:32:36
	% EndTime: 2019-10-09 23:32:36
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (4->3), mult. (16->6), div. (0->0), fcn. (16->4), ass. (0->7)
	t61 = cos(qJ(1));
	t60 = sin(qJ(1));
	t59 = cos(pkin(9));
	t58 = sin(pkin(9));
	t57 = (-t58 * t60 + t59 * t61) * qJD(1);
	t56 = (-t58 * t61 - t59 * t60) * qJD(1);
	t1 = [t56, 0, 0, 0, 0, 0; t57, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t57, 0, 0, 0, 0, 0; t56, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:32:36
	% EndTime: 2019-10-09 23:32:36
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (28->10), mult. (82->16), div. (0->0), fcn. (90->6), ass. (0->17)
	t66 = sin(qJ(5));
	t73 = qJD(5) * t66;
	t68 = cos(qJ(5));
	t72 = qJD(5) * t68;
	t64 = sin(pkin(9));
	t65 = cos(pkin(9));
	t67 = sin(qJ(1));
	t69 = cos(qJ(1));
	t62 = t69 * t64 + t67 * t65;
	t61 = -t67 * t64 + t69 * t65;
	t59 = t62 * qJD(1);
	t71 = -t59 * t66 + t61 * t72;
	t70 = -t59 * t68 - t61 * t73;
	t60 = t61 * qJD(1);
	t58 = t60 * t68 - t62 * t73;
	t57 = -t60 * t66 - t62 * t72;
	t1 = [t70, 0, 0, 0, t57, 0; t58, 0, 0, 0, t71, 0; 0, 0, 0, 0, -t73, 0; -t71, 0, 0, 0, -t58, 0; t57, 0, 0, 0, t70, 0; 0, 0, 0, 0, -t72, 0; t60, 0, 0, 0, 0, 0; t59, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:32:38
	% EndTime: 2019-10-09 23:32:38
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (108->26), mult. (317->48), div. (0->0), fcn. (353->8), ass. (0->31)
	t262 = cos(pkin(9));
	t267 = cos(qJ(1));
	t285 = sin(pkin(9));
	t286 = sin(qJ(1));
	t258 = t267 * t262 - t286 * t285;
	t263 = sin(qJ(6));
	t265 = cos(qJ(6));
	t259 = t286 * t262 + t267 * t285;
	t256 = t259 * qJD(1);
	t266 = cos(qJ(5));
	t264 = sin(qJ(5));
	t281 = qJD(5) * t264;
	t274 = t256 * t266 + t258 * t281;
	t270 = -qJD(6) * t259 + t274;
	t257 = t258 * qJD(1);
	t278 = qJD(6) * t266;
	t275 = t258 * t278 - t257;
	t287 = t270 * t263 - t275 * t265;
	t284 = t263 * t264;
	t283 = t264 * t265;
	t280 = qJD(5) * t266;
	t279 = qJD(6) * t264;
	t276 = -t259 * t278 + t256;
	t273 = t257 * t266 - t259 * t281;
	t272 = -t263 * t279 + t265 * t280;
	t271 = t263 * t280 + t265 * t279;
	t269 = qJD(6) * t258 - t273;
	t268 = -t275 * t263 - t270 * t265;
	t255 = t276 * t263 - t269 * t265;
	t254 = t269 * t263 + t276 * t265;
	t1 = [t268, 0, 0, 0, -t257 * t283 - t272 * t259, t254; t255, 0, 0, 0, -t256 * t283 + t272 * t258, -t287; 0, 0, 0, 0, -t263 * t278 - t265 * t281, -t271; t287, 0, 0, 0, t257 * t284 + t271 * t259, -t255; t254, 0, 0, 0, t256 * t284 - t271 * t258, t268; 0, 0, 0, 0, t263 * t281 - t265 * t278, -t272; -t256 * t264 + t258 * t280, 0, 0, 0, t273, 0; t257 * t264 + t259 * t280, 0, 0, 0, t274, 0; 0, 0, 0, 0, t280, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end