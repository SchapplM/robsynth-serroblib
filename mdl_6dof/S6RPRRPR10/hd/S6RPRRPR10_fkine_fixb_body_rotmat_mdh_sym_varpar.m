% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRPR10 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:49
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRRPR10_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR10_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:49:46
	% EndTime: 2020-11-04 21:49:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:49:46
	% EndTime: 2020-11-04 21:49:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t60 = cos(qJ(1));
	t59 = sin(qJ(1));
	t1 = [t60, -t59, 0, 0; t59, t60, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:49:46
	% EndTime: 2020-11-04 21:49:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t62 = cos(qJ(1));
	t61 = sin(qJ(1));
	t1 = [0, -t62, t61, t62 * pkin(1) + t61 * qJ(2) + 0; 0, -t61, -t62, t61 * pkin(1) - t62 * qJ(2) + 0; 1, 0, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:49:46
	% EndTime: 2020-11-04 21:49:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->10), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->6)
	t67 = pkin(1) + pkin(7);
	t66 = cos(qJ(1));
	t65 = cos(qJ(3));
	t64 = sin(qJ(1));
	t63 = sin(qJ(3));
	t1 = [t64 * t63, t64 * t65, t66, t64 * qJ(2) + t67 * t66 + 0; -t66 * t63, -t66 * t65, t64, -t66 * qJ(2) + t67 * t64 + 0; t65, -t63, 0, pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:49:46
	% EndTime: 2020-11-04 21:49:46
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (20->17), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->13)
	t69 = sin(qJ(4));
	t71 = sin(qJ(1));
	t79 = t71 * t69;
	t72 = cos(qJ(4));
	t78 = t71 * t72;
	t74 = cos(qJ(1));
	t77 = t74 * t69;
	t76 = t74 * t72;
	t75 = pkin(1) + pkin(7);
	t73 = cos(qJ(3));
	t70 = sin(qJ(3));
	t68 = -t70 * pkin(3) + t73 * pkin(8) - qJ(2);
	t1 = [t70 * t78 + t77, -t70 * t79 + t76, -t71 * t73, -t68 * t71 + t75 * t74 + 0; -t70 * t76 + t79, t70 * t77 + t78, t74 * t73, t68 * t74 + t75 * t71 + 0; t73 * t72, -t73 * t69, t70, t73 * pkin(3) + t70 * pkin(8) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:49:46
	% EndTime: 2020-11-04 21:49:46
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->20), mult. (36->24), div. (0->0), fcn. (49->6), ass. (0->14)
	t84 = sin(qJ(4));
	t86 = sin(qJ(1));
	t94 = t86 * t84;
	t87 = cos(qJ(4));
	t93 = t86 * t87;
	t89 = cos(qJ(1));
	t92 = t89 * t84;
	t91 = t89 * t87;
	t81 = -t87 * pkin(4) - t84 * qJ(5) - pkin(3);
	t85 = sin(qJ(3));
	t88 = cos(qJ(3));
	t90 = t88 * pkin(8) + t81 * t85 - qJ(2);
	t80 = t84 * pkin(4) - qJ(5) * t87 + pkin(1) + pkin(7);
	t1 = [t85 * t93 + t92, -t86 * t88, t85 * t94 - t91, t80 * t89 - t90 * t86 + 0; -t85 * t91 + t94, t89 * t88, -t85 * t92 - t93, t80 * t86 + t90 * t89 + 0; t88 * t87, t85, t88 * t84, t85 * pkin(8) - t81 * t88 + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:49:46
	% EndTime: 2020-11-04 21:49:46
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (48->26), mult. (56->30), div. (0->0), fcn. (79->8), ass. (0->18)
	t102 = sin(qJ(3));
	t105 = cos(qJ(3));
	t107 = pkin(8) - pkin(9);
	t101 = sin(qJ(4));
	t104 = cos(qJ(4));
	t108 = pkin(4) + pkin(5);
	t98 = t101 * qJ(5) + t108 * t104 + pkin(3);
	t115 = -t98 * t102 + t107 * t105 - qJ(2);
	t114 = cos(qJ(6));
	t113 = sin(qJ(6));
	t103 = sin(qJ(1));
	t111 = t102 * t103;
	t106 = cos(qJ(1));
	t110 = t102 * t106;
	t97 = t114 * t101 - t113 * t104;
	t96 = -t113 * t101 - t114 * t104;
	t95 = -qJ(5) * t104 + t108 * t101 + pkin(1) + pkin(7);
	t1 = [t97 * t106 - t96 * t111, t96 * t106 + t97 * t111, t103 * t105, -t115 * t103 + t95 * t106 + 0; t103 * t97 + t96 * t110, t103 * t96 - t97 * t110, -t106 * t105, t95 * t103 + t115 * t106 + 0; -t105 * t96, t105 * t97, -t102, t107 * t102 + t98 * t105 + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end