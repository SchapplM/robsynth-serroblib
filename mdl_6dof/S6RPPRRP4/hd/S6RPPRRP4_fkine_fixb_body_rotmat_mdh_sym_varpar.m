% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPRRP4 (for one body)
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
% Datum: 2020-11-04 21:28
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPPRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:28:32
	% EndTime: 2020-11-04 21:28:32
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:28:32
	% EndTime: 2020-11-04 21:28:32
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t60 = cos(qJ(1));
	t59 = sin(qJ(1));
	t1 = [t60, -t59, 0, 0; t59, t60, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:28:32
	% EndTime: 2020-11-04 21:28:32
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (6->6), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t62 = cos(qJ(1));
	t61 = sin(qJ(1));
	t1 = [t62, 0, t61, t62 * pkin(1) + t61 * qJ(2) + 0; t61, 0, -t62, t61 * pkin(1) - t62 * qJ(2) + 0; 0, 1, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:28:32
	% EndTime: 2020-11-04 21:28:32
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (12->10), mult. (12->8), div. (0->0), fcn. (20->4), ass. (0->8)
	t69 = pkin(1) + pkin(2);
	t68 = cos(qJ(1));
	t67 = sin(qJ(1));
	t66 = cos(pkin(9));
	t65 = sin(pkin(9));
	t64 = -t68 * t65 + t67 * t66;
	t63 = -t67 * t65 - t68 * t66;
	t1 = [-t63, t64, 0, t67 * qJ(2) + t69 * t68 + 0; t64, t63, 0, -t68 * qJ(2) + t69 * t67 + 0; 0, 0, -1, -qJ(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:28:32
	% EndTime: 2020-11-04 21:28:32
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (26->17), mult. (28->16), div. (0->0), fcn. (42->6), ass. (0->11)
	t80 = cos(qJ(1));
	t79 = cos(qJ(4));
	t78 = sin(qJ(1));
	t77 = sin(qJ(4));
	t76 = cos(pkin(9));
	t75 = sin(pkin(9));
	t73 = -t75 * pkin(3) + t76 * pkin(7) - qJ(2);
	t72 = t76 * pkin(3) + t75 * pkin(7) + pkin(1) + pkin(2);
	t71 = t78 * t75 + t80 * t76;
	t70 = t80 * t75 - t78 * t76;
	t1 = [t71 * t79, -t71 * t77, t70, t72 * t80 - t73 * t78 + 0; -t70 * t79, t70 * t77, t71, t72 * t78 + t73 * t80 + 0; -t77, -t79, 0, -qJ(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:28:32
	% EndTime: 2020-11-04 21:28:32
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (41->23), mult. (66->30), div. (0->0), fcn. (89->8), ass. (0->16)
	t96 = sin(qJ(1));
	t88 = sin(qJ(5));
	t91 = cos(qJ(4));
	t95 = t88 * t91;
	t90 = cos(qJ(5));
	t94 = t90 * t91;
	t89 = sin(qJ(4));
	t93 = pkin(4) * t91 + pkin(8) * t89 + pkin(3);
	t92 = cos(qJ(1));
	t87 = cos(pkin(9));
	t86 = sin(pkin(9));
	t84 = t96 * t86 + t92 * t87;
	t83 = t92 * t86 - t96 * t87;
	t82 = -t87 * pkin(7) + t93 * t86 + qJ(2);
	t81 = t86 * pkin(7) + t93 * t87 + pkin(1) + pkin(2);
	t1 = [t83 * t88 + t84 * t94, t90 * t83 - t84 * t95, t84 * t89, t81 * t92 + t82 * t96 + 0; -t83 * t94 + t84 * t88, t83 * t95 + t84 * t90, -t83 * t89, t81 * t96 - t82 * t92 + 0; -t89 * t90, t89 * t88, t91, -t89 * pkin(4) + t91 * pkin(8) + pkin(6) - qJ(3) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:28:32
	% EndTime: 2020-11-04 21:28:32
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (60->28), mult. (92->34), div. (0->0), fcn. (115->8), ass. (0->18)
	t114 = sin(qJ(1));
	t105 = sin(qJ(5));
	t108 = cos(qJ(4));
	t113 = t105 * t108;
	t107 = cos(qJ(5));
	t112 = t107 * t108;
	t101 = t107 * pkin(5) + t105 * qJ(6) + pkin(4);
	t106 = sin(qJ(4));
	t111 = pkin(8) * t106 + t101 * t108 + pkin(3);
	t110 = pkin(5) * t105 - qJ(6) * t107 + pkin(7);
	t109 = cos(qJ(1));
	t104 = cos(pkin(9));
	t103 = sin(pkin(9));
	t100 = t114 * t103 + t109 * t104;
	t99 = t109 * t103 - t114 * t104;
	t98 = t111 * t103 - t110 * t104 + qJ(2);
	t97 = t110 * t103 + t111 * t104 + pkin(1) + pkin(2);
	t1 = [t100 * t112 + t99 * t105, t100 * t106, t100 * t113 - t107 * t99, t97 * t109 + t98 * t114 + 0; t100 * t105 - t99 * t112, -t99 * t106, -t100 * t107 - t99 * t113, -t98 * t109 + t97 * t114 + 0; -t106 * t107, t108, -t106 * t105, t108 * pkin(8) - t101 * t106 + pkin(6) - qJ(3) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end