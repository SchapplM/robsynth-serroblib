% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRPP7 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:45
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRRPP7_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPP7_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:45:56
	% EndTime: 2020-11-04 21:45:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:45:56
	% EndTime: 2020-11-04 21:45:56
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t52 = cos(qJ(1));
	t51 = sin(qJ(1));
	t1 = [t52, -t51, 0, 0; t51, t52, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:45:56
	% EndTime: 2020-11-04 21:45:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t54 = cos(qJ(1));
	t53 = sin(qJ(1));
	t1 = [0, -t54, t53, t54 * pkin(1) + t53 * qJ(2) + 0; 0, -t53, -t54, t53 * pkin(1) - t54 * qJ(2) + 0; 1, 0, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:45:56
	% EndTime: 2020-11-04 21:45:56
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->10), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->6)
	t59 = pkin(1) + pkin(7);
	t58 = cos(qJ(1));
	t57 = cos(qJ(3));
	t56 = sin(qJ(1));
	t55 = sin(qJ(3));
	t1 = [t56 * t55, t56 * t57, t58, t56 * qJ(2) + t59 * t58 + 0; -t58 * t55, -t58 * t57, t56, -t58 * qJ(2) + t59 * t56 + 0; t57, -t55, 0, pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:45:56
	% EndTime: 2020-11-04 21:45:56
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (20->17), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->13)
	t61 = sin(qJ(4));
	t63 = sin(qJ(1));
	t71 = t63 * t61;
	t64 = cos(qJ(4));
	t70 = t63 * t64;
	t66 = cos(qJ(1));
	t69 = t66 * t61;
	t68 = t66 * t64;
	t67 = pkin(1) + pkin(7);
	t65 = cos(qJ(3));
	t62 = sin(qJ(3));
	t60 = -t62 * pkin(3) + t65 * pkin(8) - qJ(2);
	t1 = [t62 * t70 + t69, -t62 * t71 + t68, -t63 * t65, -t60 * t63 + t67 * t66 + 0; -t62 * t68 + t71, t62 * t69 + t70, t66 * t65, t60 * t66 + t67 * t63 + 0; t65 * t64, -t65 * t61, t62, t65 * pkin(3) + t62 * pkin(8) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:45:56
	% EndTime: 2020-11-04 21:45:56
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (29->20), mult. (36->24), div. (0->0), fcn. (49->6), ass. (0->14)
	t76 = sin(qJ(4));
	t78 = sin(qJ(1));
	t86 = t78 * t76;
	t79 = cos(qJ(4));
	t85 = t78 * t79;
	t81 = cos(qJ(1));
	t84 = t81 * t76;
	t83 = t81 * t79;
	t73 = -t79 * pkin(4) - t76 * qJ(5) - pkin(3);
	t77 = sin(qJ(3));
	t80 = cos(qJ(3));
	t82 = t80 * pkin(8) + t73 * t77 - qJ(2);
	t72 = t76 * pkin(4) - qJ(5) * t79 + pkin(1) + pkin(7);
	t1 = [t77 * t85 + t84, -t78 * t80, t77 * t86 - t83, t72 * t81 - t82 * t78 + 0; -t77 * t83 + t86, t81 * t80, -t77 * t84 - t85, t72 * t78 + t82 * t81 + 0; t80 * t79, t77, t80 * t76, t77 * pkin(8) - t73 * t80 + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:45:56
	% EndTime: 2020-11-04 21:45:56
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (38->23), mult. (36->24), div. (0->0), fcn. (49->6), ass. (0->16)
	t90 = sin(qJ(4));
	t93 = cos(qJ(4));
	t96 = pkin(4) + pkin(5);
	t88 = t90 * qJ(5) + t96 * t93 + pkin(3);
	t89 = qJ(6) - pkin(8);
	t91 = sin(qJ(3));
	t94 = cos(qJ(3));
	t103 = -t88 * t91 - t89 * t94 - qJ(2);
	t92 = sin(qJ(1));
	t101 = t92 * t90;
	t100 = t92 * t93;
	t95 = cos(qJ(1));
	t99 = t95 * t90;
	t98 = t95 * t93;
	t87 = -qJ(5) * t93 + t96 * t90 + pkin(1) + pkin(7);
	t1 = [t91 * t100 + t99, t91 * t101 - t98, t92 * t94, -t103 * t92 + t87 * t95 + 0; -t91 * t98 + t101, -t91 * t99 - t100, -t95 * t94, t103 * t95 + t87 * t92 + 0; t94 * t93, t94 * t90, -t91, t88 * t94 - t89 * t91 + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end