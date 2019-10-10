% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRP8
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:47
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPRP8_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP8_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP8_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_jacobigD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:47:49
	% EndTime: 2019-10-10 11:47:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:47:49
	% EndTime: 2019-10-10 11:47:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:47:49
	% EndTime: 2019-10-10 11:47:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, qJD(1) * cos(qJ(1)), 0, 0, 0, 0; 0, qJD(1) * sin(qJ(1)), 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:47:49
	% EndTime: 2019-10-10 11:47:49
	% DurationCPUTime: 0.09s
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
	% StartTime: 2019-10-10 11:47:49
	% EndTime: 2019-10-10 11:47:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (11->8), div. (0->0), fcn. (11->4), ass. (0->7)
	t89 = sin(qJ(1));
	t94 = qJD(1) * t89;
	t91 = cos(qJ(1));
	t93 = qJD(1) * t91;
	t92 = qJD(2) * cos(qJ(2));
	t88 = sin(qJ(2));
	t1 = [0, t93, -t88 * t94 + t91 * t92, 0, 0, 0; 0, t94, t88 * t93 + t89 * t92, 0, 0, 0; 0, 0, qJD(2) * t88, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:47:49
	% EndTime: 2019-10-10 11:47:49
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (20->8), div. (0->0), fcn. (20->4), ass. (0->10)
	t113 = sin(qJ(1));
	t119 = qJD(1) * t113;
	t115 = cos(qJ(1));
	t118 = qJD(1) * t115;
	t112 = sin(qJ(2));
	t117 = qJD(2) * t112;
	t116 = qJD(2) * cos(qJ(2));
	t111 = t112 * t118 + t113 * t116;
	t110 = t112 * t119 - t115 * t116;
	t1 = [0, t118, -t110, 0, t110, 0; 0, t119, t111, 0, -t111, 0; 0, 0, t117, 0, -t117, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:47:49
	% EndTime: 2019-10-10 11:47:49
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (20->8), div. (0->0), fcn. (20->4), ass. (0->10)
	t124 = sin(qJ(1));
	t130 = qJD(1) * t124;
	t126 = cos(qJ(1));
	t129 = qJD(1) * t126;
	t123 = sin(qJ(2));
	t128 = qJD(2) * t123;
	t127 = qJD(2) * cos(qJ(2));
	t122 = t123 * t129 + t124 * t127;
	t121 = t123 * t130 - t126 * t127;
	t1 = [0, t129, -t121, 0, t121, 0; 0, t130, t122, 0, -t122, 0; 0, 0, t128, 0, -t128, 0;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end