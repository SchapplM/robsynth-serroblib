% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPPPRR3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:29
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPPPRR3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR3_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPPRR3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t14 = qJD(1) * sin(qJ(1));
	t13 = qJD(1) * cos(qJ(1));
	t1 = [-t13, 0, 0, 0, 0, 0; -t14, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t14, 0, 0, 0, 0, 0; t13, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (4->3), mult. (16->6), div. (0->0), fcn. (16->4), ass. (0->7)
	t64 = cos(qJ(1));
	t63 = sin(qJ(1));
	t62 = cos(pkin(9));
	t61 = sin(pkin(9));
	t60 = (t61 * t64 - t62 * t63) * qJD(1);
	t59 = (t61 * t63 + t62 * t64) * qJD(1);
	t1 = [-t59, 0, 0, 0, 0, 0; t60, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t60, 0, 0, 0, 0, 0; t59, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->5), mult. (28->10), div. (0->0), fcn. (28->6), ass. (0->9)
	t44 = cos(qJ(1));
	t43 = sin(qJ(1));
	t42 = cos(pkin(9));
	t41 = cos(pkin(10));
	t40 = sin(pkin(9));
	t39 = sin(pkin(10));
	t38 = (t40 * t44 - t42 * t43) * qJD(1);
	t37 = (-t40 * t43 - t42 * t44) * qJD(1);
	t1 = [t37 * t41, 0, 0, 0, 0, 0; t38 * t41, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t37 * t39, 0, 0, 0, 0, 0; -t38 * t39, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t38, 0, 0, 0, 0, 0; t37, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (44->11), mult. (82->16), div. (0->0), fcn. (90->6), ass. (0->18)
	t74 = pkin(10) + qJ(5);
	t72 = sin(t74);
	t82 = qJD(5) * t72;
	t73 = cos(t74);
	t81 = qJD(5) * t73;
	t75 = sin(pkin(9));
	t76 = cos(pkin(9));
	t77 = sin(qJ(1));
	t78 = cos(qJ(1));
	t69 = t78 * t75 - t77 * t76;
	t70 = t77 * t75 + t78 * t76;
	t67 = t70 * qJD(1);
	t80 = -t67 * t72 + t69 * t81;
	t79 = -t67 * t73 - t69 * t82;
	t68 = t69 * qJD(1);
	t66 = t68 * t73 - t70 * t82;
	t65 = -t68 * t72 - t70 * t81;
	t1 = [t79, 0, 0, 0, t65, 0; t66, 0, 0, 0, t80, 0; 0, 0, 0, 0, t82, 0; -t80, 0, 0, 0, -t66, 0; t65, 0, 0, 0, t79, 0; 0, 0, 0, 0, t81, 0; -t68, 0, 0, 0, 0, 0; -t67, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:17
	% EndTime: 2019-10-09 23:29:17
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (162->26), mult. (317->51), div. (0->0), fcn. (353->8), ass. (0->35)
	t314 = sin(pkin(9));
	t315 = cos(pkin(9));
	t316 = sin(qJ(1));
	t317 = cos(qJ(1));
	t284 = t317 * t314 - t316 * t315;
	t292 = sin(qJ(6));
	t293 = cos(qJ(6));
	t283 = -t316 * t314 - t317 * t315;
	t281 = t283 * qJD(1);
	t291 = pkin(10) + qJ(5);
	t290 = cos(t291);
	t289 = sin(t291);
	t311 = qJD(5) * t289;
	t300 = -t281 * t290 + t284 * t311;
	t296 = -qJD(6) * t283 + t300;
	t282 = t284 * qJD(1);
	t307 = qJD(6) * t290;
	t301 = t284 * t307 + t282;
	t318 = t296 * t292 - t301 * t293;
	t313 = t289 * t292;
	t312 = t289 * t293;
	t310 = qJD(5) * t290;
	t309 = qJD(5) * t292;
	t308 = qJD(5) * t293;
	t306 = qJD(6) * t292;
	t305 = qJD(6) * t293;
	t302 = t283 * t307 + t281;
	t299 = t282 * t290 + t283 * t311;
	t298 = -t289 * t306 + t290 * t308;
	t297 = t289 * t305 + t290 * t309;
	t295 = qJD(6) * t284 + t299;
	t294 = -t301 * t292 - t296 * t293;
	t280 = t302 * t292 + t295 * t293;
	t279 = -t295 * t292 + t302 * t293;
	t1 = [t294, 0, 0, 0, -t282 * t312 + t298 * t283, t279; t280, 0, 0, 0, t281 * t312 + t298 * t284, -t318; 0, 0, 0, 0, t289 * t308 + t290 * t306, t297; t318, 0, 0, 0, t282 * t313 - t297 * t283, -t280; t279, 0, 0, 0, -t281 * t313 - t297 * t284, t294; 0, 0, 0, 0, -t289 * t309 + t290 * t305, t298; t281 * t289 + t284 * t310, 0, 0, 0, t299, 0; t282 * t289 - t283 * t310, 0, 0, 0, t300, 0; 0, 0, 0, 0, -t310, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end