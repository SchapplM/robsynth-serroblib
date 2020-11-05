% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRRP3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:51
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:51:38
	% EndTime: 2020-11-04 21:51:38
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:51:38
	% EndTime: 2020-11-04 21:51:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t64 = cos(qJ(1));
	t63 = sin(qJ(1));
	t1 = [t64, -t63, 0, 0; t63, t64, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:51:38
	% EndTime: 2020-11-04 21:51:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t67 = qJ(1) + pkin(10);
	t66 = cos(t67);
	t65 = sin(t67);
	t1 = [t66, -t65, 0, cos(qJ(1)) * pkin(1) + 0; t65, t66, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:51:38
	% EndTime: 2020-11-04 21:51:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t72 = cos(qJ(3));
	t71 = sin(qJ(3));
	t70 = qJ(1) + pkin(10);
	t69 = cos(t70);
	t68 = sin(t70);
	t1 = [t69 * t72, -t69 * t71, t68, t69 * pkin(2) + t68 * pkin(7) + cos(qJ(1)) * pkin(1) + 0; t68 * t72, -t68 * t71, -t69, t68 * pkin(2) - t69 * pkin(7) + sin(qJ(1)) * pkin(1) + 0; t71, t72, 0, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:51:38
	% EndTime: 2020-11-04 21:51:38
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (38->19), mult. (32->24), div. (0->0), fcn. (45->8), ass. (0->11)
	t76 = sin(qJ(4));
	t79 = cos(qJ(3));
	t82 = t76 * t79;
	t78 = cos(qJ(4));
	t81 = t78 * t79;
	t77 = sin(qJ(3));
	t80 = pkin(3) * t79 + pkin(8) * t77 + pkin(2);
	t75 = qJ(1) + pkin(10);
	t74 = cos(t75);
	t73 = sin(t75);
	t1 = [t73 * t76 + t74 * t81, t73 * t78 - t74 * t82, t74 * t77, cos(qJ(1)) * pkin(1) + t73 * pkin(7) + 0 + t80 * t74; t73 * t81 - t74 * t76, -t73 * t82 - t74 * t78, t73 * t77, sin(qJ(1)) * pkin(1) - t74 * pkin(7) + 0 + t80 * t73; t77 * t78, -t77 * t76, -t79, t77 * pkin(3) - t79 * pkin(8) + pkin(6) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:51:38
	% EndTime: 2020-11-04 21:51:38
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (58->23), mult. (39->26), div. (0->0), fcn. (52->10), ass. (0->15)
	t89 = qJ(4) + qJ(5);
	t86 = sin(t89);
	t92 = cos(qJ(3));
	t97 = t86 * t92;
	t87 = cos(t89);
	t96 = t87 * t92;
	t95 = pkin(4) * sin(qJ(4)) + pkin(7);
	t83 = cos(qJ(4)) * pkin(4) + pkin(3);
	t91 = sin(qJ(3));
	t93 = -pkin(9) - pkin(8);
	t94 = t83 * t92 - t91 * t93 + pkin(2);
	t88 = qJ(1) + pkin(10);
	t85 = cos(t88);
	t84 = sin(t88);
	t1 = [t84 * t86 + t85 * t96, t84 * t87 - t85 * t97, t85 * t91, cos(qJ(1)) * pkin(1) + 0 + t95 * t84 + t94 * t85; t84 * t96 - t85 * t86, -t84 * t97 - t85 * t87, t84 * t91, sin(qJ(1)) * pkin(1) + 0 - t95 * t85 + t94 * t84; t91 * t87, -t91 * t86, -t92, t91 * t83 + t92 * t93 + pkin(6) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:51:38
	% EndTime: 2020-11-04 21:51:38
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (85->28), mult. (59->32), div. (0->0), fcn. (76->10), ass. (0->19)
	t108 = qJ(4) + qJ(5);
	t105 = sin(t108);
	t111 = cos(qJ(3));
	t116 = t105 * t111;
	t106 = cos(t108);
	t115 = t106 * t111;
	t114 = pkin(4) * sin(qJ(4)) + pkin(7);
	t102 = cos(qJ(4)) * pkin(4) + pkin(3);
	t110 = sin(qJ(3));
	t112 = -pkin(9) - pkin(8);
	t113 = t102 * t111 - t110 * t112 + pkin(2);
	t107 = qJ(1) + pkin(10);
	t104 = cos(t107);
	t103 = sin(t107);
	t101 = t103 * t105 + t104 * t115;
	t100 = -t103 * t106 + t104 * t116;
	t99 = t103 * t115 - t104 * t105;
	t98 = t103 * t116 + t104 * t106;
	t1 = [t101, t104 * t110, t100, cos(qJ(1)) * pkin(1) + t101 * pkin(5) + t100 * qJ(6) + 0 + t114 * t103 + t113 * t104; t99, t103 * t110, t98, sin(qJ(1)) * pkin(1) + t99 * pkin(5) + t98 * qJ(6) + 0 - t114 * t104 + t113 * t103; t110 * t106, -t111, t110 * t105, t111 * t112 + pkin(6) + qJ(2) + 0 + (pkin(5) * t106 + qJ(6) * t105 + t102) * t110; 0, 0, 0, 1;];
	Tc_mdh = t1;
end