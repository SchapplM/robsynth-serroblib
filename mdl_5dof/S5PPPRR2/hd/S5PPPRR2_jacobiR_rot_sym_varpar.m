% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PPPRR2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PPPRR2_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPPRR2_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:59:52
	% EndTime: 2019-12-05 14:59:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:59:52
	% EndTime: 2019-12-05 14:59:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:59:52
	% EndTime: 2019-12-05 14:59:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:59:52
	% EndTime: 2019-12-05 14:59:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:59:52
	% EndTime: 2019-12-05 14:59:52
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (10->8), mult. (30->19), div. (0->0), fcn. (46->8), ass. (0->14)
	t33 = sin(pkin(8));
	t38 = sin(qJ(4));
	t42 = t33 * t38;
	t39 = cos(qJ(4));
	t41 = t33 * t39;
	t35 = cos(pkin(9));
	t36 = cos(pkin(8));
	t40 = t35 * t36;
	t37 = cos(pkin(7));
	t34 = sin(pkin(7));
	t32 = sin(pkin(9));
	t31 = t34 * t32 + t37 * t40;
	t30 = -t37 * t32 + t34 * t40;
	t1 = [0, 0, 0, -t31 * t38 + t37 * t41, 0; 0, 0, 0, -t30 * t38 + t34 * t41, 0; 0, 0, 0, -t35 * t42 - t36 * t39, 0; 0, 0, 0, -t31 * t39 - t37 * t42, 0; 0, 0, 0, -t30 * t39 - t34 * t42, 0; 0, 0, 0, -t35 * t41 + t36 * t38, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:59:52
	% EndTime: 2019-12-05 14:59:53
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (38->19), mult. (107->42), div. (0->0), fcn. (157->10), ass. (0->27)
	t91 = sin(pkin(9));
	t92 = sin(pkin(8));
	t106 = t91 * t92;
	t98 = sin(qJ(4));
	t105 = t92 * t98;
	t93 = sin(pkin(7));
	t95 = cos(pkin(8));
	t104 = t93 * t95;
	t96 = cos(pkin(7));
	t103 = t96 * t91;
	t94 = cos(pkin(9));
	t102 = t96 * t94;
	t100 = cos(qJ(4));
	t101 = t100 * t92;
	t99 = cos(qJ(5));
	t97 = sin(qJ(5));
	t90 = t94 * t101 - t95 * t98;
	t89 = -t95 * t100 - t94 * t105;
	t88 = t95 * t102 + t93 * t91;
	t87 = t95 * t103 - t93 * t94;
	t86 = t94 * t104 - t103;
	t85 = t91 * t104 + t102;
	t84 = t88 * t100 + t96 * t105;
	t83 = t96 * t101 - t88 * t98;
	t82 = t86 * t100 + t93 * t105;
	t81 = t93 * t101 - t86 * t98;
	t1 = [0, 0, 0, t83 * t99, -t84 * t97 + t87 * t99; 0, 0, 0, t81 * t99, -t82 * t97 + t85 * t99; 0, 0, 0, t89 * t99, t99 * t106 - t90 * t97; 0, 0, 0, -t83 * t97, -t84 * t99 - t87 * t97; 0, 0, 0, -t81 * t97, -t82 * t99 - t85 * t97; 0, 0, 0, -t89 * t97, -t97 * t106 - t90 * t99; 0, 0, 0, t84, 0; 0, 0, 0, t82, 0; 0, 0, 0, t90, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end