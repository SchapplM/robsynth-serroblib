% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPPR1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:33
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRPPR1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:33:04
	% EndTime: 2020-11-04 21:33:04
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:33:04
	% EndTime: 2020-11-04 21:33:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t57 = cos(qJ(1));
	t56 = sin(qJ(1));
	t1 = [t57, -t56, 0, 0; t56, t57, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:33:04
	% EndTime: 2020-11-04 21:33:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t60 = qJ(1) + pkin(9);
	t59 = cos(t60);
	t58 = sin(t60);
	t1 = [t59, -t58, 0, cos(qJ(1)) * pkin(1) + 0; t58, t59, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:33:04
	% EndTime: 2020-11-04 21:33:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t65 = cos(qJ(3));
	t64 = sin(qJ(3));
	t63 = qJ(1) + pkin(9);
	t62 = cos(t63);
	t61 = sin(t63);
	t1 = [t62 * t65, -t62 * t64, t61, t62 * pkin(2) + t61 * pkin(7) + cos(qJ(1)) * pkin(1) + 0; t61 * t65, -t61 * t64, -t62, t61 * pkin(2) - t62 * pkin(7) + sin(qJ(1)) * pkin(1) + 0; t64, t65, 0, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:33:04
	% EndTime: 2020-11-04 21:33:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (32->16), mult. (13->12), div. (0->0), fcn. (21->8), ass. (0->9)
	t73 = -qJ(4) - pkin(7);
	t72 = qJ(1) + pkin(9);
	t71 = qJ(3) + pkin(10);
	t70 = cos(t72);
	t69 = cos(t71);
	t68 = sin(t72);
	t67 = sin(t71);
	t66 = cos(qJ(3)) * pkin(3) + pkin(2);
	t1 = [t70 * t69, -t70 * t67, t68, t70 * t66 - t68 * t73 + cos(qJ(1)) * pkin(1) + 0; t68 * t69, -t68 * t67, -t70, t68 * t66 + t70 * t73 + sin(qJ(1)) * pkin(1) + 0; t67, t69, 0, sin(qJ(3)) * pkin(3) + qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:33:04
	% EndTime: 2020-11-04 21:33:04
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (58->23), mult. (35->24), div. (0->0), fcn. (48->10), ass. (0->15)
	t80 = qJ(1) + pkin(9);
	t76 = sin(t80);
	t81 = sin(pkin(11));
	t88 = t76 * t81;
	t82 = cos(pkin(11));
	t87 = t76 * t82;
	t78 = cos(t80);
	t86 = t78 * t81;
	t85 = t78 * t82;
	t79 = qJ(3) + pkin(10);
	t75 = sin(t79);
	t77 = cos(t79);
	t84 = pkin(4) * t77 + qJ(5) * t75 + cos(qJ(3)) * pkin(3) + pkin(2);
	t83 = -qJ(4) - pkin(7);
	t1 = [t77 * t85 + t88, -t77 * t86 + t87, t78 * t75, cos(qJ(1)) * pkin(1) - t76 * t83 + 0 + t84 * t78; t77 * t87 - t86, -t77 * t88 - t85, t76 * t75, sin(qJ(1)) * pkin(1) + t78 * t83 + 0 + t84 * t76; t75 * t82, -t75 * t81, -t77, t75 * pkin(4) - t77 * qJ(5) + sin(qJ(3)) * pkin(3) + qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:33:04
	% EndTime: 2020-11-04 21:33:04
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (78->27), mult. (42->27), div. (0->0), fcn. (55->12), ass. (0->17)
	t99 = qJ(1) + pkin(9);
	t93 = sin(t99);
	t98 = qJ(3) + pkin(10);
	t95 = cos(t98);
	t107 = t93 * t95;
	t97 = pkin(11) + qJ(6);
	t91 = sin(t97);
	t96 = cos(t99);
	t106 = t96 * t91;
	t94 = cos(t97);
	t105 = t96 * t94;
	t104 = sin(pkin(11)) * pkin(5) + qJ(4) + pkin(7);
	t102 = -pkin(8) - qJ(5);
	t89 = cos(pkin(11)) * pkin(5) + pkin(4);
	t92 = sin(t98);
	t103 = -t102 * t92 + t89 * t95 + cos(qJ(3)) * pkin(3) + pkin(2);
	t1 = [t95 * t105 + t93 * t91, -t95 * t106 + t93 * t94, t96 * t92, cos(qJ(1)) * pkin(1) + 0 + t104 * t93 + t103 * t96; t94 * t107 - t106, -t91 * t107 - t105, t93 * t92, sin(qJ(1)) * pkin(1) + 0 - t104 * t96 + t103 * t93; t92 * t94, -t92 * t91, -t95, t92 * t89 + t95 * t102 + sin(qJ(3)) * pkin(3) + qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end