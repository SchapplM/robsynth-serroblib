% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRRPR8_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR8_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:07:57
	% EndTime: 2019-12-29 20:07:57
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:07:57
	% EndTime: 2019-12-29 20:07:57
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0; -t8, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:07:57
	% EndTime: 2019-12-29 20:07:57
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t10 = sin(qJ(1));
	t9 = sin(qJ(2));
	t16 = t10 * t9;
	t12 = cos(qJ(1));
	t15 = t12 * t9;
	t11 = cos(qJ(2));
	t14 = t10 * t11;
	t13 = t12 * t11;
	t1 = [-t14, -t15, 0, 0, 0; t13, -t16, 0, 0, 0; 0, t11, 0, 0, 0; t16, -t13, 0, 0, 0; -t15, -t14, 0, 0, 0; 0, -t9, 0, 0, 0; t12, 0, 0, 0, 0; t10, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:07:57
	% EndTime: 2019-12-29 20:07:57
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (28->13), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t22 = qJ(2) + qJ(3);
	t20 = sin(t22);
	t23 = sin(qJ(1));
	t28 = t23 * t20;
	t21 = cos(t22);
	t27 = t23 * t21;
	t24 = cos(qJ(1));
	t26 = t24 * t20;
	t25 = t24 * t21;
	t1 = [-t27, -t26, -t26, 0, 0; t25, -t28, -t28, 0, 0; 0, t21, t21, 0, 0; t28, -t25, -t25, 0, 0; -t26, -t27, -t27, 0, 0; 0, -t20, -t20, 0, 0; t24, 0, 0, 0, 0; t23, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:07:58
	% EndTime: 2019-12-29 20:07:58
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (20->5), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t76 = cos(qJ(1));
	t75 = sin(qJ(1));
	t74 = qJ(2) + qJ(3);
	t73 = cos(t74);
	t72 = sin(t74);
	t71 = t76 * t73;
	t70 = t76 * t72;
	t69 = t75 * t73;
	t68 = t75 * t72;
	t1 = [t76, 0, 0, 0, 0; t75, 0, 0, 0, 0; 0, 0, 0, 0, 0; t69, t70, t70, 0, 0; -t71, t68, t68, 0, 0; 0, -t73, -t73, 0, 0; -t68, t71, t71, 0, 0; t70, t69, t69, 0, 0; 0, t72, t72, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:08:03
	% EndTime: 2019-12-29 20:08:03
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (44->13), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
	t104 = qJ(2) + qJ(3);
	t102 = sin(t104);
	t106 = sin(qJ(1));
	t114 = t106 * t102;
	t105 = sin(qJ(5));
	t113 = t106 * t105;
	t107 = cos(qJ(5));
	t112 = t106 * t107;
	t108 = cos(qJ(1));
	t111 = t108 * t102;
	t110 = t108 * t105;
	t109 = t108 * t107;
	t103 = cos(t104);
	t101 = t102 * t107;
	t100 = t102 * t105;
	t99 = t103 * t109;
	t98 = t103 * t110;
	t97 = t103 * t112;
	t96 = t103 * t113;
	t95 = -t102 * t113 + t109;
	t94 = t102 * t112 + t110;
	t93 = t102 * t110 + t112;
	t92 = t102 * t109 - t113;
	t1 = [t95, t98, t98, 0, t92; t93, t96, t96, 0, t94; 0, t100, t100, 0, -t103 * t107; -t94, t99, t99, 0, -t93; t92, t97, t97, 0, t95; 0, t101, t101, 0, t103 * t105; -t106 * t103, -t111, -t111, 0, 0; t108 * t103, -t114, -t114, 0, 0; 0, t103, t103, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end