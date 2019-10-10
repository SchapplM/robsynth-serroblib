% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:28
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRPR4_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR4_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0, 0; -t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->5)
	t10 = cos(qJ(1));
	t9 = sin(qJ(1));
	t8 = cos(pkin(10));
	t7 = sin(pkin(10));
	t1 = [-t9 * t8, 0, 0, 0, 0, 0; t10 * t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t9 * t7, 0, 0, 0, 0, 0; -t10 * t7, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t10, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->8), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t14 = pkin(10) + qJ(3);
	t12 = sin(t14);
	t15 = sin(qJ(1));
	t20 = t15 * t12;
	t13 = cos(t14);
	t19 = t15 * t13;
	t16 = cos(qJ(1));
	t18 = t16 * t12;
	t17 = t16 * t13;
	t1 = [-t19, 0, -t18, 0, 0, 0; t17, 0, -t20, 0, 0, 0; 0, 0, t13, 0, 0, 0; t20, 0, -t17, 0, 0, 0; -t18, 0, -t19, 0, 0, 0; 0, 0, -t12, 0, 0, 0; t16, 0, 0, 0, 0, 0; t15, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (44->13), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t23 = pkin(10) + qJ(3) + qJ(4);
	t21 = sin(t23);
	t24 = sin(qJ(1));
	t29 = t24 * t21;
	t22 = cos(t23);
	t28 = t24 * t22;
	t25 = cos(qJ(1));
	t27 = t25 * t21;
	t26 = t25 * t22;
	t1 = [-t28, 0, -t27, -t27, 0, 0; t26, 0, -t29, -t29, 0, 0; 0, 0, t22, t22, 0, 0; t29, 0, -t26, -t26, 0, 0; -t27, 0, -t28, -t28, 0, 0; 0, 0, -t21, -t21, 0, 0; t25, 0, 0, 0, 0, 0; t24, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (59->12), mult. (38->18), div. (0->0), fcn. (66->6), ass. (0->20)
	t88 = pkin(10) + qJ(3) + qJ(4);
	t87 = cos(t88);
	t89 = sin(pkin(11));
	t99 = t87 * t89;
	t91 = sin(qJ(1));
	t98 = t91 * t89;
	t90 = cos(pkin(11));
	t97 = t91 * t90;
	t92 = cos(qJ(1));
	t96 = t92 * t89;
	t95 = t92 * t90;
	t86 = sin(t88);
	t94 = t86 * t97;
	t93 = t86 * t95;
	t85 = t92 * t87;
	t84 = t91 * t87;
	t83 = t87 * t90;
	t82 = t86 * t96;
	t81 = t86 * t98;
	t1 = [-t87 * t97 + t96, 0, -t93, -t93, 0, 0; t87 * t95 + t98, 0, -t94, -t94, 0, 0; 0, 0, t83, t83, 0, 0; t87 * t98 + t95, 0, t82, t82, 0, 0; -t87 * t96 + t97, 0, t81, t81, 0, 0; 0, 0, -t99, -t99, 0, 0; -t91 * t86, 0, t85, t85, 0, 0; t92 * t86, 0, t84, t84, 0, 0; 0, 0, t86, t86, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (107->17), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->25)
	t103 = pkin(10) + qJ(3) + qJ(4);
	t100 = cos(t103);
	t104 = pkin(11) + qJ(6);
	t101 = sin(t104);
	t113 = t100 * t101;
	t105 = sin(qJ(1));
	t112 = t105 * t101;
	t102 = cos(t104);
	t111 = t105 * t102;
	t106 = cos(qJ(1));
	t110 = t106 * t101;
	t109 = t106 * t102;
	t99 = sin(t103);
	t108 = t99 * t111;
	t107 = t99 * t109;
	t98 = t106 * t100;
	t97 = t105 * t100;
	t96 = t100 * t102;
	t95 = t99 * t110;
	t94 = t99 * t112;
	t93 = t100 * t109 + t112;
	t92 = -t100 * t110 + t111;
	t91 = -t100 * t111 + t110;
	t90 = t100 * t112 + t109;
	t1 = [t91, 0, -t107, -t107, 0, t92; t93, 0, -t108, -t108, 0, -t90; 0, 0, t96, t96, 0, -t99 * t101; t90, 0, t95, t95, 0, -t93; t92, 0, t94, t94, 0, t91; 0, 0, -t113, -t113, 0, -t99 * t102; -t105 * t99, 0, t98, t98, 0, 0; t106 * t99, 0, t97, t97, 0, 0; 0, 0, t99, t99, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end