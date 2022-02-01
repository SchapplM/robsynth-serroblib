% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPRR6
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
%   Siehe auch: S5RRPRR6_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRPRR6_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR6_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
JR_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:18:05
	% EndTime: 2022-01-20 11:18:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:18:05
	% EndTime: 2022-01-20 11:18:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0; -t8, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:18:05
	% EndTime: 2022-01-20 11:18:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (14->7), mult. (0->0), div. (0->0), fcn. (8->2), ass. (0->4)
	t16 = qJ(1) + qJ(2);
	t15 = cos(t16);
	t14 = sin(t16);
	t1 = [-t14, -t14, 0, 0, 0; t15, t15, 0, 0, 0; 0, 0, 0, 0, 0; -t15, -t15, 0, 0, 0; -t14, -t14, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:18:05
	% EndTime: 2022-01-20 11:18:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (16->5), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t23 = qJ(1) + qJ(2);
	t21 = sin(t23);
	t25 = cos(pkin(9));
	t27 = t21 * t25;
	t22 = cos(t23);
	t24 = sin(pkin(9));
	t26 = t22 * t24;
	t20 = t22 * t25;
	t19 = t21 * t24;
	t1 = [-t27, -t27, 0, 0, 0; t20, t20, 0, 0, 0; 0, 0, 0, 0, 0; t19, t19, 0, 0, 0; -t26, -t26, 0, 0, 0; 0, 0, 0, 0, 0; t22, t22, 0, 0, 0; t21, t21, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:18:05
	% EndTime: 2022-01-20 11:18:05
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (44->11), mult. (42->14), div. (0->0), fcn. (72->6), ass. (0->16)
	t69 = qJ(1) + qJ(2);
	t67 = sin(t69);
	t70 = sin(pkin(9));
	t76 = t67 * t70;
	t71 = cos(pkin(9));
	t72 = sin(qJ(4));
	t75 = t71 * t72;
	t73 = cos(qJ(4));
	t74 = t71 * t73;
	t68 = cos(t69);
	t66 = t68 * t70;
	t65 = t67 * t72 + t68 * t74;
	t64 = t67 * t73 - t68 * t75;
	t63 = -t67 * t74 + t68 * t72;
	t62 = t67 * t75 + t68 * t73;
	t1 = [t63, t63, 0, t64, 0; t65, t65, 0, -t62, 0; 0, 0, 0, -t70 * t72, 0; t62, t62, 0, -t65, 0; t64, t64, 0, t63, 0; 0, 0, 0, -t70 * t73, 0; -t76, -t76, 0, 0, 0; t66, t66, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:18:06
	% EndTime: 2022-01-20 11:18:06
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (94->16), mult. (56->14), div. (0->0), fcn. (96->6), ass. (0->19)
	t92 = qJ(1) + qJ(2);
	t88 = sin(t92);
	t93 = sin(pkin(9));
	t99 = t88 * t93;
	t94 = cos(pkin(9));
	t98 = t88 * t94;
	t90 = cos(t92);
	t97 = t90 * t94;
	t91 = qJ(4) + qJ(5);
	t87 = sin(t91);
	t96 = t93 * t87;
	t89 = cos(t91);
	t95 = t93 * t89;
	t86 = t90 * t93;
	t85 = t88 * t87 + t89 * t97;
	t84 = -t87 * t97 + t88 * t89;
	t83 = t90 * t87 - t89 * t98;
	t82 = t87 * t98 + t90 * t89;
	t1 = [t83, t83, 0, t84, t84; t85, t85, 0, -t82, -t82; 0, 0, 0, -t96, -t96; t82, t82, 0, -t85, -t85; t84, t84, 0, t83, t83; 0, 0, 0, -t95, -t95; -t99, -t99, 0, 0, 0; t86, t86, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
end