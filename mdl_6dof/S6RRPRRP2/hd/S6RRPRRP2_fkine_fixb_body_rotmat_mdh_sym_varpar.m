% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRP2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:13
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPRRP2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:13:28
	% EndTime: 2020-11-04 22:13:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:13:28
	% EndTime: 2020-11-04 22:13:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t62 = cos(qJ(1));
	t61 = sin(qJ(1));
	t1 = [t62, -t61, 0, 0; t61, t62, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:13:28
	% EndTime: 2020-11-04 22:13:28
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t66 = cos(qJ(1));
	t65 = cos(qJ(2));
	t64 = sin(qJ(1));
	t63 = sin(qJ(2));
	t1 = [t66 * t65, -t66 * t63, t64, t66 * pkin(1) + t64 * pkin(7) + 0; t64 * t65, -t64 * t63, -t66, t64 * pkin(1) - t66 * pkin(7) + 0; t63, t65, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:13:28
	% EndTime: 2020-11-04 22:13:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t73 = cos(qJ(1));
	t72 = sin(qJ(1));
	t71 = -qJ(3) - pkin(7);
	t70 = qJ(2) + pkin(10);
	t69 = cos(t70);
	t68 = sin(t70);
	t67 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t73 * t69, -t73 * t68, t72, t73 * t67 - t71 * t72 + 0; t72 * t69, -t72 * t68, -t73, t72 * t67 + t73 * t71 + 0; t68, t69, 0, sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:13:28
	% EndTime: 2020-11-04 22:13:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->15), mult. (14->12), div. (0->0), fcn. (22->8), ass. (0->9)
	t79 = qJ(2) + pkin(10);
	t81 = cos(qJ(1));
	t80 = sin(qJ(1));
	t78 = -pkin(8) - qJ(3) - pkin(7);
	t77 = qJ(4) + t79;
	t76 = cos(t77);
	t75 = sin(t77);
	t74 = pkin(3) * cos(t79) + cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t81 * t76, -t81 * t75, t80, t81 * t74 - t80 * t78 + 0; t80 * t76, -t80 * t75, -t81, t80 * t74 + t81 * t78 + 0; t75, t76, 0, pkin(3) * sin(t79) + sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:13:28
	% EndTime: 2020-11-04 22:13:28
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (60->22), mult. (36->24), div. (0->0), fcn. (49->10), ass. (0->15)
	t88 = sin(qJ(5));
	t89 = sin(qJ(1));
	t96 = t89 * t88;
	t90 = cos(qJ(5));
	t95 = t89 * t90;
	t91 = cos(qJ(1));
	t94 = t91 * t88;
	t93 = t91 * t90;
	t87 = qJ(2) + pkin(10);
	t85 = qJ(4) + t87;
	t83 = sin(t85);
	t84 = cos(t85);
	t92 = pkin(4) * t84 + pkin(9) * t83 + pkin(3) * cos(t87) + cos(qJ(2)) * pkin(2) + pkin(1);
	t86 = -pkin(8) - qJ(3) - pkin(7);
	t1 = [t84 * t93 + t96, -t84 * t94 + t95, t91 * t83, -t89 * t86 + t92 * t91 + 0; t84 * t95 - t94, -t84 * t96 - t93, t89 * t83, t91 * t86 + t92 * t89 + 0; t83 * t90, -t83 * t88, -t84, t83 * pkin(4) - t84 * pkin(9) + pkin(3) * sin(t87) + sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:13:28
	% EndTime: 2020-11-04 22:13:28
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (81->27), mult. (56->30), div. (0->0), fcn. (73->10), ass. (0->19)
	t107 = sin(qJ(5));
	t108 = sin(qJ(1));
	t115 = t108 * t107;
	t109 = cos(qJ(5));
	t114 = t108 * t109;
	t110 = cos(qJ(1));
	t113 = t110 * t107;
	t112 = t110 * t109;
	t106 = qJ(2) + pkin(10);
	t104 = qJ(4) + t106;
	t102 = sin(t104);
	t103 = cos(t104);
	t111 = pkin(4) * t103 + pkin(9) * t102 + pkin(3) * cos(t106) + cos(qJ(2)) * pkin(2) + pkin(1);
	t105 = -pkin(8) - qJ(3) - pkin(7);
	t100 = t103 * t112 + t115;
	t99 = t103 * t113 - t114;
	t98 = t103 * t114 - t113;
	t97 = t103 * t115 + t112;
	t1 = [t100, t110 * t102, t99, t100 * pkin(5) + t99 * qJ(6) - t108 * t105 + t111 * t110 + 0; t98, t108 * t102, t97, t98 * pkin(5) + t97 * qJ(6) + t110 * t105 + t111 * t108 + 0; t102 * t109, -t103, t102 * t107, pkin(3) * sin(t106) + sin(qJ(2)) * pkin(2) - t103 * pkin(9) + pkin(6) + 0 + (pkin(5) * t109 + qJ(6) * t107 + pkin(4)) * t102; 0, 0, 0, 1;];
	Tc_mdh = t1;
end