% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRRP1 (for one body)
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
% Datum: 2020-11-04 21:50
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRRRP1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:50:57
	% EndTime: 2020-11-04 21:50:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:50:57
	% EndTime: 2020-11-04 21:50:57
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t60 = cos(qJ(1));
	t59 = sin(qJ(1));
	t1 = [t60, -t59, 0, 0; t59, t60, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:50:57
	% EndTime: 2020-11-04 21:50:57
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t63 = qJ(1) + pkin(10);
	t62 = cos(t63);
	t61 = sin(t63);
	t1 = [t62, -t61, 0, cos(qJ(1)) * pkin(1) + 0; t61, t62, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:50:57
	% EndTime: 2020-11-04 21:50:57
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t68 = cos(qJ(3));
	t67 = sin(qJ(3));
	t66 = qJ(1) + pkin(10);
	t65 = cos(t66);
	t64 = sin(t66);
	t1 = [t65 * t68, -t65 * t67, t64, t65 * pkin(2) + t64 * pkin(7) + cos(qJ(1)) * pkin(1) + 0; t64 * t68, -t64 * t67, -t65, t64 * pkin(2) - t65 * pkin(7) + sin(qJ(1)) * pkin(1) + 0; t67, t68, 0, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:50:57
	% EndTime: 2020-11-04 21:50:57
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (32->16), mult. (13->12), div. (0->0), fcn. (21->8), ass. (0->9)
	t76 = -pkin(8) - pkin(7);
	t75 = qJ(3) + qJ(4);
	t74 = qJ(1) + pkin(10);
	t73 = cos(t75);
	t72 = sin(t75);
	t71 = cos(t74);
	t70 = sin(t74);
	t69 = cos(qJ(3)) * pkin(3) + pkin(2);
	t1 = [t71 * t73, -t71 * t72, t70, t71 * t69 - t70 * t76 + cos(qJ(1)) * pkin(1) + 0; t70 * t73, -t70 * t72, -t71, t70 * t69 + t71 * t76 + sin(qJ(1)) * pkin(1) + 0; t72, t73, 0, sin(qJ(3)) * pkin(3) + qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:50:57
	% EndTime: 2020-11-04 21:50:57
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (58->23), mult. (35->26), div. (0->0), fcn. (48->10), ass. (0->13)
	t83 = qJ(3) + qJ(4);
	t81 = cos(t83);
	t84 = sin(qJ(5));
	t89 = t81 * t84;
	t85 = cos(qJ(5));
	t88 = t81 * t85;
	t80 = sin(t83);
	t87 = pkin(4) * t81 + pkin(9) * t80 + cos(qJ(3)) * pkin(3) + pkin(2);
	t86 = -pkin(8) - pkin(7);
	t82 = qJ(1) + pkin(10);
	t79 = cos(t82);
	t78 = sin(t82);
	t1 = [t78 * t84 + t79 * t88, t78 * t85 - t79 * t89, t79 * t80, cos(qJ(1)) * pkin(1) - t78 * t86 + 0 + t87 * t79; t78 * t88 - t79 * t84, -t78 * t89 - t79 * t85, t78 * t80, sin(qJ(1)) * pkin(1) + t79 * t86 + 0 + t87 * t78; t80 * t85, -t80 * t84, -t81, t80 * pkin(4) - t81 * pkin(9) + sin(qJ(3)) * pkin(3) + qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:50:57
	% EndTime: 2020-11-04 21:50:57
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (81->28), mult. (55->30), div. (0->0), fcn. (72->10), ass. (0->19)
	t101 = sin(qJ(5));
	t99 = qJ(1) + pkin(10);
	t95 = sin(t99);
	t108 = t95 * t101;
	t102 = cos(qJ(5));
	t107 = t95 * t102;
	t96 = cos(t99);
	t106 = t96 * t101;
	t105 = t96 * t102;
	t100 = qJ(3) + qJ(4);
	t97 = sin(t100);
	t98 = cos(t100);
	t104 = pkin(4) * t98 + pkin(9) * t97 + cos(qJ(3)) * pkin(3) + pkin(2);
	t103 = -pkin(8) - pkin(7);
	t93 = t98 * t105 + t108;
	t92 = t98 * t106 - t107;
	t91 = t98 * t107 - t106;
	t90 = t98 * t108 + t105;
	t1 = [t93, t96 * t97, t92, cos(qJ(1)) * pkin(1) + t93 * pkin(5) + t92 * qJ(6) - t95 * t103 + 0 + t104 * t96; t91, t95 * t97, t90, sin(qJ(1)) * pkin(1) + t91 * pkin(5) + t90 * qJ(6) + t96 * t103 + 0 + t104 * t95; t97 * t102, -t98, t97 * t101, sin(qJ(3)) * pkin(3) - t98 * pkin(9) + pkin(6) + qJ(2) + 0 + (pkin(5) * t102 + qJ(6) * t101 + pkin(4)) * t97; 0, 0, 0, 1;];
	Tc_mdh = t1;
end