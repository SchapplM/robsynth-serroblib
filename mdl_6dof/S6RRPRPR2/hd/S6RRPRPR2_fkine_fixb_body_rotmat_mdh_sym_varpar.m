% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRPR2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:08
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPRPR2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:08:18
	% EndTime: 2020-11-04 22:08:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:08:18
	% EndTime: 2020-11-04 22:08:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t57 = cos(qJ(1));
	t56 = sin(qJ(1));
	t1 = [t57, -t56, 0, 0; t56, t57, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:08:18
	% EndTime: 2020-11-04 22:08:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t61 = cos(qJ(1));
	t60 = cos(qJ(2));
	t59 = sin(qJ(1));
	t58 = sin(qJ(2));
	t1 = [t61 * t60, -t61 * t58, t59, t61 * pkin(1) + t59 * pkin(7) + 0; t59 * t60, -t59 * t58, -t61, t59 * pkin(1) - t61 * pkin(7) + 0; t58, t60, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:08:18
	% EndTime: 2020-11-04 22:08:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t68 = cos(qJ(1));
	t67 = sin(qJ(1));
	t66 = -qJ(3) - pkin(7);
	t65 = qJ(2) + pkin(10);
	t64 = cos(t65);
	t63 = sin(t65);
	t62 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t68 * t64, -t68 * t63, t67, t68 * t62 - t66 * t67 + 0; t67 * t64, -t67 * t63, -t68, t67 * t62 + t68 * t66 + 0; t63, t64, 0, sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:08:18
	% EndTime: 2020-11-04 22:08:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->15), mult. (14->12), div. (0->0), fcn. (22->8), ass. (0->9)
	t74 = qJ(2) + pkin(10);
	t76 = cos(qJ(1));
	t75 = sin(qJ(1));
	t73 = -pkin(8) - qJ(3) - pkin(7);
	t72 = qJ(4) + t74;
	t71 = cos(t72);
	t70 = sin(t72);
	t69 = pkin(3) * cos(t74) + cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t76 * t71, -t76 * t70, t75, t76 * t69 - t75 * t73 + 0; t75 * t71, -t75 * t70, -t76, t75 * t69 + t76 * t73 + 0; t70, t71, 0, pkin(3) * sin(t74) + sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:08:18
	% EndTime: 2020-11-04 22:08:18
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (53->21), mult. (24->16), div. (0->0), fcn. (32->8), ass. (0->9)
	t82 = qJ(2) + pkin(10);
	t80 = qJ(4) + t82;
	t78 = sin(t80);
	t79 = cos(t80);
	t85 = pkin(4) * t79 + qJ(5) * t78 + pkin(3) * cos(t82) + cos(qJ(2)) * pkin(2) + pkin(1);
	t84 = cos(qJ(1));
	t83 = sin(qJ(1));
	t81 = -pkin(8) - qJ(3) - pkin(7);
	t1 = [t83, -t84 * t79, t84 * t78, -t83 * t81 + t85 * t84 + 0; -t84, -t83 * t79, t83 * t78, t84 * t81 + t85 * t83 + 0; 0, -t78, -t79, t78 * pkin(4) - t79 * qJ(5) + pkin(3) * sin(t82) + sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:08:18
	% EndTime: 2020-11-04 22:08:18
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (71->23), mult. (43->24), div. (0->0), fcn. (56->10), ass. (0->16)
	t102 = pkin(4) + pkin(9);
	t101 = pkin(5) + pkin(8) + qJ(3) + pkin(7);
	t92 = sin(qJ(6));
	t93 = sin(qJ(1));
	t100 = t93 * t92;
	t94 = cos(qJ(6));
	t99 = t93 * t94;
	t95 = cos(qJ(1));
	t98 = t95 * t92;
	t97 = t95 * t94;
	t91 = qJ(2) + pkin(10);
	t89 = qJ(4) + t91;
	t87 = sin(t89);
	t88 = cos(t89);
	t96 = qJ(5) * t87 + t102 * t88 + pkin(3) * cos(t91) + cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t87 * t98 + t99, t87 * t97 - t100, t95 * t88, t101 * t93 + t96 * t95 + 0; t87 * t100 - t97, t87 * t99 + t98, t93 * t88, -t101 * t95 + t96 * t93 + 0; -t88 * t92, -t88 * t94, t87, pkin(3) * sin(t91) + sin(qJ(2)) * pkin(2) - t88 * qJ(5) + pkin(6) + 0 + t102 * t87; 0, 0, 0, 1;];
	Tc_mdh = t1;
end