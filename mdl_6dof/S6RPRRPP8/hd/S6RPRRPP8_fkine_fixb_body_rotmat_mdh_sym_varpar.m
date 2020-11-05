% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRPP8 (for one body)
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
% Datum: 2020-11-04 21:46
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRRPP8_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPP8_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:46:15
	% EndTime: 2020-11-04 21:46:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:46:15
	% EndTime: 2020-11-04 21:46:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t53 = cos(qJ(1));
	t52 = sin(qJ(1));
	t1 = [t53, -t52, 0, 0; t52, t53, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:46:16
	% EndTime: 2020-11-04 21:46:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t55 = cos(qJ(1));
	t54 = sin(qJ(1));
	t1 = [0, -t55, t54, t55 * pkin(1) + t54 * qJ(2) + 0; 0, -t54, -t55, t54 * pkin(1) - t55 * qJ(2) + 0; 1, 0, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:46:16
	% EndTime: 2020-11-04 21:46:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->10), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->6)
	t60 = pkin(1) + pkin(7);
	t59 = cos(qJ(1));
	t58 = cos(qJ(3));
	t57 = sin(qJ(1));
	t56 = sin(qJ(3));
	t1 = [t57 * t56, t57 * t58, t59, t57 * qJ(2) + t60 * t59 + 0; -t59 * t56, -t59 * t58, t57, -t59 * qJ(2) + t60 * t57 + 0; t58, -t56, 0, pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:46:16
	% EndTime: 2020-11-04 21:46:16
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (20->17), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->13)
	t62 = sin(qJ(4));
	t64 = sin(qJ(1));
	t72 = t64 * t62;
	t65 = cos(qJ(4));
	t71 = t64 * t65;
	t67 = cos(qJ(1));
	t70 = t67 * t62;
	t69 = t67 * t65;
	t68 = pkin(1) + pkin(7);
	t66 = cos(qJ(3));
	t63 = sin(qJ(3));
	t61 = -t63 * pkin(3) + t66 * pkin(8) - qJ(2);
	t1 = [t63 * t71 + t70, -t63 * t72 + t69, -t64 * t66, -t61 * t64 + t68 * t67 + 0; -t63 * t69 + t72, t63 * t70 + t71, t67 * t66, t61 * t67 + t68 * t64 + 0; t66 * t65, -t66 * t62, t63, t66 * pkin(3) + t63 * pkin(8) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:46:16
	% EndTime: 2020-11-04 21:46:16
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (30->21), mult. (36->24), div. (0->0), fcn. (49->6), ass. (0->14)
	t77 = sin(qJ(4));
	t79 = sin(qJ(1));
	t87 = t79 * t77;
	t80 = cos(qJ(4));
	t86 = t79 * t80;
	t82 = cos(qJ(1));
	t85 = t82 * t77;
	t84 = t82 * t80;
	t74 = -t80 * pkin(4) - t77 * qJ(5) - pkin(3);
	t78 = sin(qJ(3));
	t81 = cos(qJ(3));
	t83 = t81 * pkin(8) + t74 * t78 - qJ(2);
	t73 = t77 * pkin(4) - qJ(5) * t80 + pkin(1) + pkin(7);
	t1 = [-t79 * t81, -t78 * t86 - t85, t78 * t87 - t84, t73 * t82 - t83 * t79 + 0; t82 * t81, t78 * t84 - t87, -t78 * t85 - t86, t73 * t79 + t83 * t82 + 0; t78, -t81 * t80, t81 * t77, t78 * pkin(8) - t74 * t81 + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:46:16
	% EndTime: 2020-11-04 21:46:16
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (37->22), mult. (36->24), div. (0->0), fcn. (49->6), ass. (0->16)
	t90 = qJ(6) + pkin(4);
	t91 = sin(qJ(4));
	t94 = cos(qJ(4));
	t89 = qJ(5) * t91 + t90 * t94 + pkin(3);
	t92 = sin(qJ(3));
	t95 = cos(qJ(3));
	t97 = pkin(5) + pkin(8);
	t104 = -t89 * t92 + t95 * t97 - qJ(2);
	t93 = sin(qJ(1));
	t102 = t93 * t91;
	t101 = t93 * t94;
	t96 = cos(qJ(1));
	t100 = t96 * t91;
	t99 = t96 * t94;
	t88 = -qJ(5) * t94 + t90 * t91 + pkin(1) + pkin(7);
	t1 = [-t93 * t95, t102 * t92 - t99, t101 * t92 + t100, -t104 * t93 + t88 * t96 + 0; t96 * t95, -t100 * t92 - t101, -t92 * t99 + t102, t104 * t96 + t88 * t93 + 0; t92, t95 * t91, t95 * t94, t89 * t95 + t92 * t97 + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end