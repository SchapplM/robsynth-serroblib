% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:07
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRPR3_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR3_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR3_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_jacobigD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:07:53
	% EndTime: 2019-10-10 10:07:53
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:07:53
	% EndTime: 2019-10-10 10:07:53
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:07:53
	% EndTime: 2019-10-10 10:07:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, qJD(1) * cos(qJ(1)), 0, 0, 0, 0; 0, qJD(1) * sin(qJ(1)), 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:07:53
	% EndTime: 2019-10-10 10:07:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, qJD(1) * cos(qJ(1)), 0, 0, 0, 0; 0, qJD(1) * sin(qJ(1)), 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:07:53
	% EndTime: 2019-10-10 10:07:53
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (7->3), mult. (11->8), div. (0->0), fcn. (11->4), ass. (0->8)
	t85 = sin(qJ(1));
	t89 = qJD(1) * t85;
	t86 = cos(qJ(1));
	t88 = qJD(1) * t86;
	t84 = qJ(2) + pkin(10);
	t87 = qJD(2) * cos(t84);
	t82 = sin(t84);
	t1 = [0, t88, 0, -t82 * t89 + t86 * t87, 0, 0; 0, t89, 0, t82 * t88 + t85 * t87, 0, 0; 0, 0, 0, qJD(2) * t82, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:07:53
	% EndTime: 2019-10-10 10:07:53
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (7->3), mult. (11->8), div. (0->0), fcn. (11->4), ass. (0->8)
	t94 = sin(qJ(1));
	t98 = qJD(1) * t94;
	t95 = cos(qJ(1));
	t97 = qJD(1) * t95;
	t93 = qJ(2) + pkin(10);
	t96 = qJD(2) * cos(t93);
	t91 = sin(t93);
	t1 = [0, t97, 0, -t91 * t98 + t95 * t96, 0, 0; 0, t98, 0, t91 * t97 + t94 * t96, 0, 0; 0, 0, 0, qJD(2) * t91, 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:07:53
	% EndTime: 2019-10-10 10:07:53
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (14->3), mult. (20->8), div. (0->0), fcn. (20->4), ass. (0->11)
	t111 = sin(qJ(1));
	t115 = qJD(1) * t111;
	t112 = cos(qJ(1));
	t114 = qJD(1) * t112;
	t110 = qJ(2) + pkin(10);
	t113 = qJD(2) * cos(t110);
	t108 = sin(t110);
	t107 = qJD(2) * t108;
	t106 = t108 * t114 + t111 * t113;
	t105 = -t108 * t115 + t112 * t113;
	t1 = [0, t114, 0, t105, 0, t105; 0, t115, 0, t106, 0, t106; 0, 0, 0, t107, 0, t107;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end