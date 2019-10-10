% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:51
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRPRR4_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR4_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR4_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_jacobigD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:51:27
	% EndTime: 2019-10-10 00:51:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:51:27
	% EndTime: 2019-10-10 00:51:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:51:27
	% EndTime: 2019-10-10 00:51:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:51:27
	% EndTime: 2019-10-10 00:51:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->1), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->2)
	t20 = qJ(1) + pkin(10);
	t1 = [0, 0, qJD(1) * cos(t20), 0, 0, 0; 0, 0, qJD(1) * sin(t20), 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:51:27
	% EndTime: 2019-10-10 00:51:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->1), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->2)
	t51 = qJ(1) + pkin(10);
	t1 = [0, 0, qJD(1) * cos(t51), 0, 0, 0; 0, 0, qJD(1) * sin(t51), 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:51:27
	% EndTime: 2019-10-10 00:51:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->3), mult. (11->9), div. (0->0), fcn. (11->4), ass. (0->7)
	t79 = cos(qJ(3));
	t81 = qJD(1) * t79;
	t80 = qJD(3) * sin(qJ(3));
	t77 = qJ(1) + pkin(10);
	t76 = cos(t77);
	t75 = sin(t77);
	t1 = [0, 0, qJD(1) * t76, 0, -t75 * t81 - t76 * t80, 0; 0, 0, qJD(1) * t75, 0, -t75 * t80 + t76 * t81, 0; 0, 0, 0, 0, qJD(3) * t79, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:51:27
	% EndTime: 2019-10-10 00:51:27
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (14->3), mult. (20->9), div. (0->0), fcn. (20->4), ass. (0->10)
	t104 = cos(qJ(3));
	t106 = qJD(1) * t104;
	t105 = qJD(3) * sin(qJ(3));
	t102 = qJ(1) + pkin(10);
	t101 = qJD(3) * t104;
	t100 = cos(t102);
	t99 = sin(t102);
	t98 = t100 * t106 - t99 * t105;
	t97 = -t100 * t105 - t99 * t106;
	t1 = [0, 0, qJD(1) * t100, 0, t97, t97; 0, 0, qJD(1) * t99, 0, t98, t98; 0, 0, 0, 0, t101, t101;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end