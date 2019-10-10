% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR7
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:46
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPPRR7_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR7_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_jacobig_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t43 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t43, 0, 0, 0, 0; 0, -cos(qJ(1)) * t43, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t54 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t54, 0, 0, 0, 0; 0, -cos(qJ(1)) * t54, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t49 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t49, 0, 0, 0, 0; 0, -cos(qJ(1)) * t49, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->3), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t69 = cos(pkin(6));
	t70 = sin(qJ(2));
	t74 = t69 * t70;
	t73 = cos(qJ(1));
	t72 = cos(qJ(2));
	t71 = sin(qJ(1));
	t68 = sin(pkin(6));
	t1 = [0, t71 * t68, 0, 0, -t71 * t74 + t73 * t72, 0; 0, -t73 * t68, 0, 0, t71 * t72 + t73 * t74, 0; 1, t69, 0, 0, t68 * t70, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (8->8), mult. (24->19), div. (0->0), fcn. (40->8), ass. (0->14)
	t103 = cos(qJ(5));
	t98 = sin(pkin(6));
	t110 = t98 * t103;
	t101 = sin(qJ(2));
	t102 = sin(qJ(1));
	t109 = t102 * t101;
	t104 = cos(qJ(2));
	t108 = t102 * t104;
	t105 = cos(qJ(1));
	t107 = t105 * t101;
	t106 = t105 * t104;
	t100 = sin(qJ(5));
	t99 = cos(pkin(6));
	t1 = [0, t102 * t98, 0, 0, -t99 * t109 + t106, (t99 * t108 + t107) * t100 + t102 * t110; 0, -t105 * t98, 0, 0, t99 * t107 + t108, (-t99 * t106 + t109) * t100 - t105 * t110; 1, t99, 0, 0, t98 * t101, -t98 * t104 * t100 + t99 * t103;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end