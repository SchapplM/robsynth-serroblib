% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRRP2 (for one body)
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

function Tc_mdh = S6RPRRRP2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:51:18
	% EndTime: 2020-11-04 21:51:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:51:18
	% EndTime: 2020-11-04 21:51:18
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t59 = cos(qJ(1));
	t58 = sin(qJ(1));
	t1 = [t59, -t58, 0, 0; t58, t59, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:51:18
	% EndTime: 2020-11-04 21:51:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t62 = qJ(1) + pkin(10);
	t61 = cos(t62);
	t60 = sin(t62);
	t1 = [t61, -t60, 0, cos(qJ(1)) * pkin(1) + 0; t60, t61, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:51:18
	% EndTime: 2020-11-04 21:51:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t67 = cos(qJ(3));
	t66 = sin(qJ(3));
	t65 = qJ(1) + pkin(10);
	t64 = cos(t65);
	t63 = sin(t65);
	t1 = [t64 * t67, -t64 * t66, t63, t64 * pkin(2) + t63 * pkin(7) + cos(qJ(1)) * pkin(1) + 0; t63 * t67, -t63 * t66, -t64, t63 * pkin(2) - t64 * pkin(7) + sin(qJ(1)) * pkin(1) + 0; t66, t67, 0, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:51:18
	% EndTime: 2020-11-04 21:51:18
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (38->19), mult. (32->24), div. (0->0), fcn. (45->8), ass. (0->11)
	t71 = sin(qJ(4));
	t74 = cos(qJ(3));
	t77 = t71 * t74;
	t73 = cos(qJ(4));
	t76 = t73 * t74;
	t72 = sin(qJ(3));
	t75 = pkin(3) * t74 + pkin(8) * t72 + pkin(2);
	t70 = qJ(1) + pkin(10);
	t69 = cos(t70);
	t68 = sin(t70);
	t1 = [t68 * t71 + t69 * t76, t68 * t73 - t69 * t77, t69 * t72, cos(qJ(1)) * pkin(1) + t68 * pkin(7) + 0 + t75 * t69; t68 * t76 - t69 * t71, -t68 * t77 - t69 * t73, t68 * t72, sin(qJ(1)) * pkin(1) - t69 * pkin(7) + 0 + t75 * t68; t72 * t73, -t72 * t71, -t74, t72 * pkin(3) - t74 * pkin(8) + pkin(6) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:51:18
	% EndTime: 2020-11-04 21:51:18
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (58->23), mult. (39->26), div. (0->0), fcn. (52->10), ass. (0->15)
	t84 = qJ(4) + qJ(5);
	t81 = sin(t84);
	t87 = cos(qJ(3));
	t92 = t81 * t87;
	t82 = cos(t84);
	t91 = t82 * t87;
	t90 = pkin(4) * sin(qJ(4)) + pkin(7);
	t78 = cos(qJ(4)) * pkin(4) + pkin(3);
	t86 = sin(qJ(3));
	t88 = -pkin(9) - pkin(8);
	t89 = t78 * t87 - t86 * t88 + pkin(2);
	t83 = qJ(1) + pkin(10);
	t80 = cos(t83);
	t79 = sin(t83);
	t1 = [t79 * t81 + t80 * t91, t79 * t82 - t80 * t92, t80 * t86, cos(qJ(1)) * pkin(1) + 0 + t90 * t79 + t89 * t80; t79 * t91 - t80 * t81, -t79 * t92 - t80 * t82, t79 * t86, sin(qJ(1)) * pkin(1) + 0 - t90 * t80 + t89 * t79; t86 * t82, -t86 * t81, -t87, t86 * t78 + t87 * t88 + pkin(6) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:51:18
	% EndTime: 2020-11-04 21:51:18
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (71->25), mult. (44->28), div. (0->0), fcn. (57->10), ass. (0->15)
	t101 = qJ(4) + qJ(5);
	t97 = sin(t101);
	t107 = pkin(7) + pkin(5) * t97 + pkin(4) * sin(qJ(4));
	t103 = cos(qJ(3));
	t106 = t103 * t97;
	t98 = cos(t101);
	t105 = t103 * t98;
	t102 = sin(qJ(3));
	t93 = pkin(5) * t98 + cos(qJ(4)) * pkin(4) + pkin(3);
	t99 = -qJ(6) - pkin(9) - pkin(8);
	t104 = -t102 * t99 + t103 * t93 + pkin(2);
	t100 = qJ(1) + pkin(10);
	t96 = cos(t100);
	t95 = sin(t100);
	t1 = [t96 * t105 + t95 * t97, -t96 * t106 + t95 * t98, t96 * t102, cos(qJ(1)) * pkin(1) + 0 + t107 * t95 + t104 * t96; t95 * t105 - t96 * t97, -t95 * t106 - t96 * t98, t95 * t102, sin(qJ(1)) * pkin(1) + 0 - t107 * t96 + t104 * t95; t102 * t98, -t102 * t97, -t103, t102 * t93 + t103 * t99 + pkin(6) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end