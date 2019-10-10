% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:02
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPRR6_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR6_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR6_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_jacobigD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:02:21
	% EndTime: 2019-10-10 12:02:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:02:21
	% EndTime: 2019-10-10 12:02:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:02:21
	% EndTime: 2019-10-10 12:02:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, qJD(1) * cos(qJ(1)), 0, 0, 0, 0; 0, qJD(1) * sin(qJ(1)), 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:02:21
	% EndTime: 2019-10-10 12:02:21
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (2->2), mult. (11->8), div. (0->0), fcn. (11->4), ass. (0->7)
	t73 = sin(qJ(1));
	t78 = qJD(1) * t73;
	t75 = cos(qJ(1));
	t77 = qJD(1) * t75;
	t76 = qJD(2) * cos(qJ(2));
	t72 = sin(qJ(2));
	t1 = [0, t77, -t72 * t78 + t75 * t76, 0, 0, 0; 0, t78, t72 * t77 + t73 * t76, 0, 0, 0; 0, 0, qJD(2) * t72, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:02:21
	% EndTime: 2019-10-10 12:02:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (11->8), div. (0->0), fcn. (11->4), ass. (0->7)
	t80 = sin(qJ(1));
	t85 = qJD(1) * t80;
	t82 = cos(qJ(1));
	t84 = qJD(1) * t82;
	t83 = qJD(2) * cos(qJ(2));
	t79 = sin(qJ(2));
	t1 = [0, t84, -t79 * t85 + t82 * t83, 0, 0, 0; 0, t85, t79 * t84 + t80 * t83, 0, 0, 0; 0, 0, qJD(2) * t79, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:02:21
	% EndTime: 2019-10-10 12:02:21
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->2), mult. (20->8), div. (0->0), fcn. (20->4), ass. (0->10)
	t99 = sin(qJ(1));
	t104 = qJD(1) * t99;
	t101 = cos(qJ(1));
	t103 = qJD(1) * t101;
	t102 = qJD(2) * cos(qJ(2));
	t98 = sin(qJ(2));
	t97 = qJD(2) * t98;
	t96 = t99 * t102 + t98 * t103;
	t95 = t101 * t102 - t98 * t104;
	t1 = [0, t103, t95, 0, t95, 0; 0, t104, t96, 0, t96, 0; 0, 0, t97, 0, t97, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:02:21
	% EndTime: 2019-10-10 12:02:21
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (6->2), mult. (29->8), div. (0->0), fcn. (29->4), ass. (0->10)
	t105 = sin(qJ(1));
	t110 = qJD(1) * t105;
	t107 = cos(qJ(1));
	t109 = qJD(1) * t107;
	t108 = qJD(2) * cos(qJ(2));
	t104 = sin(qJ(2));
	t103 = qJD(2) * t104;
	t102 = t104 * t109 + t105 * t108;
	t101 = -t104 * t110 + t107 * t108;
	t1 = [0, t109, t101, 0, t101, t101; 0, t110, t102, 0, t102, t102; 0, 0, t103, 0, t103, t103;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end