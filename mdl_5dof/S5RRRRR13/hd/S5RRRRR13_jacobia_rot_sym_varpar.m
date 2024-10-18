% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRRR13
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5RRRRR13_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RRRRR13_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR13_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
Ja_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (18->0), mult. (12->0), div. (10->0), fcn. (12->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (45->0), mult. (18->0), div. (15->0), fcn. (18->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 1, 1, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (383->16), mult. (294->35), div. (65->9), fcn. (455->9), ass. (0->27)
	t52 = cos(pkin(5));
	t48 = qJ(1) + qJ(2) + qJ(3);
	t47 = cos(t48);
	t51 = sin(pkin(5));
	t59 = t47 * t51;
	t45 = atan2(t59, t52);
	t42 = sin(t45);
	t43 = cos(t45);
	t38 = t42 * t59 + t43 * t52;
	t46 = sin(t48);
	t61 = 0.1e1 / t38 ^ 2 * t46 ^ 2;
	t49 = t51 ^ 2;
	t44 = 0.1e1 / (0.1e1 + t47 ^ 2 * t49 / t52 ^ 2);
	t60 = t44 / t52;
	t53 = sin(qJ(4));
	t58 = t52 * t53;
	t54 = cos(qJ(4));
	t57 = t52 * t54;
	t41 = -t46 * t58 + t47 * t54;
	t39 = 0.1e1 / t41 ^ 2;
	t40 = t46 * t57 + t47 * t53;
	t56 = t40 ^ 2 * t39 + 0.1e1;
	t55 = t46 * t51 * t60;
	t36 = 0.1e1 / t56;
	t34 = ((-t46 * t53 + t47 * t57) / t41 - (-t46 * t54 - t47 * t58) * t40 * t39) * t36;
	t33 = (0.1e1 / t38 * t59 - (-t43 * t47 * t49 * t60 + (t44 - 0.1e1) * t51 * t42) * t51 * t61) / (t49 * t61 + 0.1e1);
	t1 = [-t55, -t55, -t55, 0, 0; t33, t33, t33, 0, 0; t34, t34, t34, t56 * t36, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (711->22), mult. (378->38), div. (72->9), fcn. (510->13), ass. (0->32)
	t83 = cos(pkin(5));
	t78 = qJ(1) + qJ(2) + qJ(3);
	t73 = cos(t78);
	t82 = sin(pkin(5));
	t85 = t73 * t82;
	t70 = atan2(t85, t83);
	t65 = sin(t70);
	t66 = cos(t70);
	t60 = t65 * t85 + t66 * t83;
	t72 = sin(t78);
	t88 = 0.1e1 / t60 ^ 2 * t72 ^ 2;
	t81 = qJ(4) + qJ(5);
	t74 = pkin(5) + t81;
	t75 = pkin(5) - t81;
	t68 = sin(t74) / 0.2e1 - sin(t75) / 0.2e1;
	t77 = cos(t81);
	t64 = -t72 * t68 + t73 * t77;
	t62 = 0.1e1 / t64 ^ 2;
	t69 = cos(t75) / 0.2e1 + cos(t74) / 0.2e1;
	t76 = sin(t81);
	t63 = t72 * t69 + t73 * t76;
	t87 = t63 ^ 2 * t62;
	t79 = t82 ^ 2;
	t67 = 0.1e1 / (0.1e1 + t73 ^ 2 * t79 / t83 ^ 2);
	t86 = t67 / t83;
	t84 = t72 * t82 * t86;
	t61 = 0.1e1 / t64;
	t57 = 0.1e1 / (0.1e1 + t87);
	t56 = (t64 * t61 + t87) * t57;
	t55 = ((t73 * t69 - t72 * t76) * t61 - (-t73 * t68 - t72 * t77) * t63 * t62) * t57;
	t54 = (0.1e1 / t60 * t85 - (-t66 * t73 * t79 * t86 + (t67 - 0.1e1) * t82 * t65) * t82 * t88) / (t79 * t88 + 0.1e1);
	t1 = [-t84, -t84, -t84, 0, 0; t54, t54, t54, 0, 0; t55, t55, t55, t56, t56;];
	Ja_rot = t1;
end