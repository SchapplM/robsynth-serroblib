% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPRRP8 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:29
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPPRRP8_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRP8_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:29:55
	% EndTime: 2020-11-04 21:29:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:29:55
	% EndTime: 2020-11-04 21:29:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t57 = cos(qJ(1));
	t56 = sin(qJ(1));
	t1 = [t57, -t56, 0, 0; t56, t57, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:29:55
	% EndTime: 2020-11-04 21:29:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t59 = cos(qJ(1));
	t58 = sin(qJ(1));
	t1 = [0, -t59, t58, t59 * pkin(1) + t58 * qJ(2) + 0; 0, -t58, -t59, t58 * pkin(1) - t59 * qJ(2) + 0; 1, 0, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:29:55
	% EndTime: 2020-11-04 21:29:55
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (11->10), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->6)
	t64 = cos(qJ(1));
	t63 = sin(qJ(1));
	t62 = pkin(1) + qJ(3);
	t61 = cos(pkin(9));
	t60 = sin(pkin(9));
	t1 = [t63 * t60, t63 * t61, t64, qJ(2) * t63 + t62 * t64 + 0; -t64 * t60, -t64 * t61, t63, -qJ(2) * t64 + t62 * t63 + 0; t61, -t60, 0, pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:29:55
	% EndTime: 2020-11-04 21:29:55
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->13), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t71 = cos(qJ(1));
	t70 = sin(qJ(1));
	t69 = pkin(9) + qJ(4);
	t68 = pkin(1) + pkin(7) + qJ(3);
	t67 = cos(t69);
	t66 = sin(t69);
	t65 = sin(pkin(9)) * pkin(3) + qJ(2);
	t1 = [t70 * t66, t70 * t67, t71, t65 * t70 + t68 * t71 + 0; -t71 * t66, -t71 * t67, t70, -t65 * t71 + t68 * t70 + 0; t67, -t66, 0, cos(pkin(9)) * pkin(3) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:29:55
	% EndTime: 2020-11-04 21:29:55
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (40->20), mult. (33->22), div. (0->0), fcn. (46->8), ass. (0->14)
	t77 = sin(qJ(5));
	t78 = sin(qJ(1));
	t85 = t78 * t77;
	t79 = cos(qJ(5));
	t84 = t78 * t79;
	t80 = cos(qJ(1));
	t83 = t80 * t77;
	t82 = t80 * t79;
	t76 = pkin(9) + qJ(4);
	t73 = sin(t76);
	t74 = cos(t76);
	t81 = pkin(4) * t73 - pkin(8) * t74 + sin(pkin(9)) * pkin(3) + qJ(2);
	t75 = pkin(1) + pkin(7) + qJ(3);
	t1 = [t73 * t84 + t83, -t73 * t85 + t82, -t78 * t74, t75 * t80 + t81 * t78 + 0; -t73 * t82 + t85, t73 * t83 + t84, t80 * t74, t75 * t78 - t81 * t80 + 0; t74 * t79, -t74 * t77, t73, t74 * pkin(4) + t73 * pkin(8) + cos(pkin(9)) * pkin(3) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:29:55
	% EndTime: 2020-11-04 21:29:55
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (55->26), mult. (53->28), div. (0->0), fcn. (70->8), ass. (0->18)
	t95 = sin(qJ(5));
	t96 = sin(qJ(1));
	t103 = t96 * t95;
	t97 = cos(qJ(5));
	t102 = t96 * t97;
	t98 = cos(qJ(1));
	t101 = t98 * t95;
	t100 = t98 * t97;
	t94 = pkin(9) + qJ(4);
	t91 = sin(t94);
	t92 = cos(t94);
	t99 = pkin(4) * t91 - pkin(8) * t92 + sin(pkin(9)) * pkin(3) + qJ(2);
	t93 = pkin(1) + pkin(7) + qJ(3);
	t89 = -t91 * t100 + t103;
	t88 = t91 * t101 + t102;
	t87 = t91 * t102 + t101;
	t86 = t91 * t103 - t100;
	t1 = [t87, -t96 * t92, t86, t87 * pkin(5) + t86 * qJ(6) + t93 * t98 + t99 * t96 + 0; t89, t98 * t92, -t88, t89 * pkin(5) - t88 * qJ(6) + t93 * t96 - t99 * t98 + 0; t92 * t97, t91, t92 * t95, cos(pkin(9)) * pkin(3) + t91 * pkin(8) + pkin(2) + pkin(6) + 0 + (pkin(5) * t97 + qJ(6) * t95 + pkin(4)) * t92; 0, 0, 0, 1;];
	Tc_mdh = t1;
end