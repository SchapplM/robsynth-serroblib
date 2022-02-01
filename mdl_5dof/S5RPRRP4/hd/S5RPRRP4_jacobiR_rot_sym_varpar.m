% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRRP4
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
%   Siehe auch: S5RPRRP4_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPRRP4_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP4_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
JR_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0; -t8, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
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
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->9), mult. (28->12), div. (0->0), fcn. (48->6), ass. (0->15)
	t51 = sin(qJ(3));
	t52 = sin(qJ(1));
	t58 = t52 * t51;
	t53 = cos(qJ(3));
	t57 = t52 * t53;
	t54 = cos(qJ(1));
	t56 = t54 * t51;
	t55 = t54 * t53;
	t50 = cos(pkin(8));
	t49 = sin(pkin(8));
	t48 = t50 * t55 + t58;
	t47 = -t50 * t56 + t57;
	t46 = -t50 * t57 + t56;
	t45 = t50 * t58 + t55;
	t1 = [t46, 0, t47, 0, 0; t48, 0, -t45, 0, 0; 0, 0, -t49 * t51, 0, 0; t45, 0, -t48, 0, 0; t47, 0, t46, 0, 0; 0, 0, -t49 * t53, 0, 0; -t52 * t49, 0, 0, 0, 0; t54 * t49, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (45->14), mult. (42->12), div. (0->0), fcn. (72->6), ass. (0->18)
	t67 = qJ(3) + qJ(4);
	t65 = sin(t67);
	t68 = sin(pkin(8));
	t77 = t68 * t65;
	t66 = cos(t67);
	t76 = t68 * t66;
	t70 = sin(qJ(1));
	t75 = t70 * t65;
	t74 = t70 * t66;
	t71 = cos(qJ(1));
	t73 = t71 * t65;
	t72 = t71 * t66;
	t69 = cos(pkin(8));
	t64 = t69 * t72 + t75;
	t63 = -t69 * t73 + t74;
	t62 = -t69 * t74 + t73;
	t61 = t69 * t75 + t72;
	t1 = [t62, 0, t63, t63, 0; t64, 0, -t61, -t61, 0; 0, 0, -t77, -t77, 0; t61, 0, -t64, -t64, 0; t63, 0, t62, t62, 0; 0, 0, -t76, -t76, 0; -t70 * t68, 0, 0, 0, 0; t71 * t68, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (45->14), mult. (42->12), div. (0->0), fcn. (72->6), ass. (0->18)
	t71 = qJ(3) + qJ(4);
	t69 = sin(t71);
	t72 = sin(pkin(8));
	t81 = t72 * t69;
	t70 = cos(t71);
	t80 = t72 * t70;
	t74 = sin(qJ(1));
	t79 = t74 * t69;
	t78 = t74 * t70;
	t75 = cos(qJ(1));
	t77 = t75 * t69;
	t76 = t75 * t70;
	t73 = cos(pkin(8));
	t68 = t73 * t76 + t79;
	t67 = -t73 * t77 + t78;
	t66 = -t73 * t78 + t77;
	t65 = t73 * t79 + t76;
	t1 = [t66, 0, t67, t67, 0; t68, 0, -t65, -t65, 0; 0, 0, -t81, -t81, 0; t65, 0, -t68, -t68, 0; t67, 0, t66, t66, 0; 0, 0, -t80, -t80, 0; -t74 * t72, 0, 0, 0, 0; t75 * t72, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
end