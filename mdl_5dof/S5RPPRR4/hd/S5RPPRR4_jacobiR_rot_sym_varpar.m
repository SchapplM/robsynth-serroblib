% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPPRR4
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
%   Siehe auch: S5RPPRR4_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPPRR4_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR4_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
JR_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:41
	% EndTime: 2022-01-23 09:17:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:41
	% EndTime: 2022-01-23 09:17:41
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
	% StartTime: 2022-01-23 09:17:41
	% EndTime: 2022-01-23 09:17:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->5)
	t10 = cos(qJ(1));
	t9 = sin(qJ(1));
	t8 = cos(pkin(8));
	t7 = sin(pkin(8));
	t1 = [-t9 * t8, 0, 0, 0, 0; t10 * t8, 0, 0, 0, 0; 0, 0, 0, 0, 0; t9 * t7, 0, 0, 0, 0; -t10 * t7, 0, 0, 0, 0; 0, 0, 0, 0, 0; t10, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:41
	% EndTime: 2022-01-23 09:17:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (5->5), mult. (14->10), div. (0->0), fcn. (24->6), ass. (0->11)
	t38 = sin(pkin(9));
	t42 = sin(qJ(1));
	t47 = t42 * t38;
	t40 = cos(pkin(9));
	t46 = t42 * t40;
	t43 = cos(qJ(1));
	t45 = t43 * t38;
	t44 = t43 * t40;
	t41 = cos(pkin(8));
	t39 = sin(pkin(8));
	t1 = [-t41 * t46 + t45, 0, 0, 0, 0; t41 * t44 + t47, 0, 0, 0, 0; 0, 0, 0, 0, 0; t41 * t47 + t44, 0, 0, 0, 0; -t41 * t45 + t46, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t39 * t42, 0, 0, 0, 0; t39 * t43, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:41
	% EndTime: 2022-01-23 09:17:41
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (29->10), mult. (28->12), div. (0->0), fcn. (48->6), ass. (0->16)
	t52 = pkin(9) + qJ(4);
	t50 = sin(t52);
	t55 = sin(qJ(1));
	t60 = t55 * t50;
	t51 = cos(t52);
	t59 = t55 * t51;
	t56 = cos(qJ(1));
	t58 = t56 * t50;
	t57 = t56 * t51;
	t54 = cos(pkin(8));
	t53 = sin(pkin(8));
	t49 = t54 * t57 + t60;
	t48 = -t54 * t58 + t59;
	t47 = -t54 * t59 + t58;
	t46 = t54 * t60 + t57;
	t1 = [t47, 0, 0, t48, 0; t49, 0, 0, -t46, 0; 0, 0, 0, -t53 * t50, 0; t46, 0, 0, -t49, 0; t48, 0, 0, t47, 0; 0, 0, 0, -t53 * t51, 0; -t53 * t55, 0, 0, 0, 0; t53 * t56, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:41
	% EndTime: 2022-01-23 09:17:41
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (73->14), mult. (42->12), div. (0->0), fcn. (72->6), ass. (0->18)
	t72 = pkin(9) + qJ(4) + qJ(5);
	t70 = sin(t72);
	t73 = sin(pkin(8));
	t82 = t73 * t70;
	t71 = cos(t72);
	t81 = t73 * t71;
	t75 = sin(qJ(1));
	t80 = t75 * t70;
	t79 = t75 * t71;
	t76 = cos(qJ(1));
	t78 = t76 * t70;
	t77 = t76 * t71;
	t74 = cos(pkin(8));
	t69 = t74 * t77 + t80;
	t68 = -t74 * t78 + t79;
	t67 = -t74 * t79 + t78;
	t66 = t74 * t80 + t77;
	t1 = [t67, 0, 0, t68, t68; t69, 0, 0, -t66, -t66; 0, 0, 0, -t82, -t82; t66, 0, 0, -t69, -t69; t68, 0, 0, t67, t67; 0, 0, 0, -t81, -t81; -t75 * t73, 0, 0, 0, 0; t76 * t73, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
end