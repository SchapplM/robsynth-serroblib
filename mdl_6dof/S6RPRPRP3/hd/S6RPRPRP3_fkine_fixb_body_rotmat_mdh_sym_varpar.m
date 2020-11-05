% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPRP3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:36
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRPRP3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:36:28
	% EndTime: 2020-11-04 21:36:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:36:28
	% EndTime: 2020-11-04 21:36:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t62 = cos(qJ(1));
	t61 = sin(qJ(1));
	t1 = [t62, -t61, 0, 0; t61, t62, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:36:28
	% EndTime: 2020-11-04 21:36:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t65 = qJ(1) + pkin(9);
	t64 = cos(t65);
	t63 = sin(t65);
	t1 = [t64, -t63, 0, cos(qJ(1)) * pkin(1) + 0; t63, t64, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:36:28
	% EndTime: 2020-11-04 21:36:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t70 = cos(qJ(3));
	t69 = sin(qJ(3));
	t68 = qJ(1) + pkin(9);
	t67 = cos(t68);
	t66 = sin(t68);
	t1 = [t67 * t70, -t67 * t69, t66, t67 * pkin(2) + t66 * pkin(7) + cos(qJ(1)) * pkin(1) + 0; t66 * t70, -t66 * t69, -t67, t66 * pkin(2) - t67 * pkin(7) + sin(qJ(1)) * pkin(1) + 0; t69, t70, 0, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:36:28
	% EndTime: 2020-11-04 21:36:28
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (38->19), mult. (32->24), div. (0->0), fcn. (45->8), ass. (0->11)
	t74 = sin(pkin(10));
	t77 = cos(qJ(3));
	t80 = t74 * t77;
	t75 = cos(pkin(10));
	t79 = t75 * t77;
	t76 = sin(qJ(3));
	t78 = pkin(3) * t77 + qJ(4) * t76 + pkin(2);
	t73 = qJ(1) + pkin(9);
	t72 = cos(t73);
	t71 = sin(t73);
	t1 = [t71 * t74 + t72 * t79, t71 * t75 - t72 * t80, t72 * t76, cos(qJ(1)) * pkin(1) + t71 * pkin(7) + 0 + t78 * t72; t71 * t79 - t72 * t74, -t71 * t80 - t72 * t75, t71 * t76, sin(qJ(1)) * pkin(1) - t72 * pkin(7) + 0 + t78 * t71; t76 * t75, -t76 * t74, -t77, t76 * pkin(3) - t77 * qJ(4) + pkin(6) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:36:28
	% EndTime: 2020-11-04 21:36:28
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (58->23), mult. (39->26), div. (0->0), fcn. (52->10), ass. (0->15)
	t87 = qJ(1) + pkin(9);
	t83 = sin(t87);
	t91 = cos(qJ(3));
	t95 = t83 * t91;
	t85 = cos(t87);
	t94 = t85 * t91;
	t93 = sin(pkin(10)) * pkin(4) + pkin(7);
	t81 = cos(pkin(10)) * pkin(4) + pkin(3);
	t89 = -pkin(8) - qJ(4);
	t90 = sin(qJ(3));
	t92 = t81 * t91 - t89 * t90 + pkin(2);
	t86 = pkin(10) + qJ(5);
	t84 = cos(t86);
	t82 = sin(t86);
	t1 = [t83 * t82 + t84 * t94, -t82 * t94 + t83 * t84, t85 * t90, cos(qJ(1)) * pkin(1) + 0 + t93 * t83 + t92 * t85; -t85 * t82 + t84 * t95, -t82 * t95 - t85 * t84, t83 * t90, sin(qJ(1)) * pkin(1) + 0 - t93 * t85 + t92 * t83; t90 * t84, -t90 * t82, -t91, t90 * t81 + t91 * t89 + pkin(6) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:36:28
	% EndTime: 2020-11-04 21:36:28
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (85->28), mult. (59->32), div. (0->0), fcn. (76->10), ass. (0->19)
	t106 = qJ(1) + pkin(9);
	t102 = sin(t106);
	t110 = cos(qJ(3));
	t114 = t102 * t110;
	t104 = cos(t106);
	t113 = t104 * t110;
	t112 = sin(pkin(10)) * pkin(4) + pkin(7);
	t100 = cos(pkin(10)) * pkin(4) + pkin(3);
	t108 = -pkin(8) - qJ(4);
	t109 = sin(qJ(3));
	t111 = t100 * t110 - t108 * t109 + pkin(2);
	t105 = pkin(10) + qJ(5);
	t103 = cos(t105);
	t101 = sin(t105);
	t99 = t102 * t101 + t103 * t113;
	t98 = t101 * t113 - t102 * t103;
	t97 = -t104 * t101 + t103 * t114;
	t96 = t101 * t114 + t104 * t103;
	t1 = [t99, t104 * t109, t98, cos(qJ(1)) * pkin(1) + t99 * pkin(5) + t98 * qJ(6) + 0 + t112 * t102 + t111 * t104; t97, t102 * t109, t96, sin(qJ(1)) * pkin(1) + t97 * pkin(5) + t96 * qJ(6) + 0 - t112 * t104 + t111 * t102; t109 * t103, -t110, t109 * t101, t110 * t108 + pkin(6) + qJ(2) + 0 + (pkin(5) * t103 + qJ(6) * t101 + t100) * t109; 0, 0, 0, 1;];
	Tc_mdh = t1;
end