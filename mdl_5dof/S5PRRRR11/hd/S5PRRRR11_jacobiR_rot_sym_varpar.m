% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRRRR11
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
%   Siehe auch: S5PRRRR11_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PRRRR11_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR11_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
JR_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:36
	% EndTime: 2024-09-27 21:46:36
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	% Symbolic code from jacobiR_rot_1_floatb_twist_matlab.m not found
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:36
	% EndTime: 2024-09-27 21:46:36
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (7->4), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->4)
	t12 = pkin(10) + qJ(2);
	t11 = cos(t12);
	t10 = sin(t12);
	t1 = [0, -t10, 0, 0, 0; 0, t11, 0, 0, 0; 0, 0, 0, 0, 0; 0, -t11, 0, 0, 0; 0, -t10, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:36
	% EndTime: 2024-09-27 21:46:36
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (27->8), mult. (28->14), div. (0->0), fcn. (48->6), ass. (0->14)
	t53 = cos(pkin(5));
	t54 = sin(qJ(3));
	t57 = t53 * t54;
	t55 = cos(qJ(3));
	t56 = t53 * t55;
	t52 = sin(pkin(5));
	t51 = pkin(10) + qJ(2);
	t50 = cos(t51);
	t49 = sin(t51);
	t48 = -t49 * t57 + t50 * t55;
	t47 = -t49 * t56 - t50 * t54;
	t46 = -t49 * t55 - t50 * t57;
	t45 = t49 * t54 - t50 * t56;
	t1 = [0, t46, t47, 0, 0; 0, t48, -t45, 0, 0; 0, 0, t52 * t55, 0, 0; 0, t45, -t48, 0, 0; 0, t47, t46, 0, 0; 0, 0, -t52 * t54, 0, 0; 0, t50 * t52, 0, 0, 0; 0, t49 * t52, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:36
	% EndTime: 2024-09-27 21:46:36
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (130->15), mult. (58->16), div. (0->0), fcn. (72->9), ass. (0->22)
	t94 = qJ(3) + qJ(4);
	t90 = pkin(5) + t94;
	t86 = cos(t90) / 0.2e1;
	t97 = pkin(5) - t94;
	t87 = cos(t97);
	t99 = t87 / 0.2e1 + t86;
	t98 = sin(t90) / 0.2e1;
	t96 = sin(t97);
	t82 = t98 - t96 / 0.2e1;
	t93 = pkin(10) + qJ(2);
	t88 = sin(t93);
	t89 = cos(t93);
	t92 = cos(t94);
	t77 = -t89 * t82 - t88 * t92;
	t79 = t88 * t82 - t89 * t92;
	t95 = sin(pkin(5));
	t91 = sin(t94);
	t83 = t86 - t87 / 0.2e1;
	t81 = t98 + t96 / 0.2e1;
	t78 = -t88 * t99 - t89 * t91;
	t76 = t88 * t91 - t89 * t99;
	t1 = [0, t77, t78, t78, 0; 0, -t79, -t76, -t76, 0; 0, 0, t81, t81, 0; 0, t76, t79, t79, 0; 0, t78, t77, t77, 0; 0, 0, t83, t83, 0; 0, t89 * t95, 0, 0, 0; 0, t88 * t95, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:37
	% EndTime: 2024-09-27 21:46:37
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (236->18), mult. (78->16), div. (0->0), fcn. (96->9), ass. (0->22)
	t101 = qJ(3) + qJ(4) + qJ(5);
	t96 = pkin(5) + t101;
	t110 = cos(t96) / 0.2e1;
	t109 = sin(t96) / 0.2e1;
	t108 = pkin(5) - t101;
	t102 = pkin(10) + qJ(2);
	t100 = cos(t102);
	t105 = sin(t108);
	t92 = t109 - t105 / 0.2e1;
	t98 = cos(t101);
	t99 = sin(t102);
	t107 = t100 * t98 - t99 * t92;
	t87 = -t100 * t92 - t99 * t98;
	t106 = cos(t108);
	t104 = t106 / 0.2e1 + t110;
	t103 = sin(pkin(5));
	t97 = sin(t101);
	t93 = t110 - t106 / 0.2e1;
	t91 = t109 + t105 / 0.2e1;
	t88 = -t100 * t97 - t99 * t104;
	t86 = -t100 * t104 + t99 * t97;
	t1 = [0, t87, t88, t88, t88; 0, t107, -t86, -t86, -t86; 0, 0, t91, t91, t91; 0, t86, -t107, -t107, -t107; 0, t88, t87, t87, t87; 0, 0, t93, t93, t93; 0, t100 * t103, 0, 0, 0; 0, t99 * t103, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
end