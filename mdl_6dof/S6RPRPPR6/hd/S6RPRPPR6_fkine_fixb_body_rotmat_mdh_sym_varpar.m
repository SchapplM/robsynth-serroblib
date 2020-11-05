% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPPR6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:34
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRPPR6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:34:46
	% EndTime: 2020-11-04 21:34:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:34:46
	% EndTime: 2020-11-04 21:34:46
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t60 = cos(qJ(1));
	t59 = sin(qJ(1));
	t1 = [t60, -t59, 0, 0; t59, t60, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:34:46
	% EndTime: 2020-11-04 21:34:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t62 = cos(qJ(1));
	t61 = sin(qJ(1));
	t1 = [0, -t62, t61, t62 * pkin(1) + t61 * qJ(2) + 0; 0, -t61, -t62, t61 * pkin(1) - t62 * qJ(2) + 0; 1, 0, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:34:46
	% EndTime: 2020-11-04 21:34:46
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
	% StartTime: 2020-11-04 21:34:46
	% EndTime: 2020-11-04 21:34:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->13), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t74 = cos(qJ(1));
	t73 = sin(qJ(1));
	t72 = qJ(3) + pkin(9);
	t71 = pkin(1) + pkin(7) + qJ(4);
	t70 = cos(t72);
	t69 = sin(t72);
	t68 = sin(qJ(3)) * pkin(3) + qJ(2);
	t1 = [t73 * t69, t73 * t70, t74, t68 * t73 + t71 * t74 + 0; -t74 * t69, -t74 * t70, t73, -t68 * t74 + t71 * t73 + 0; t70, -t69, 0, cos(qJ(3)) * pkin(3) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:34:46
	% EndTime: 2020-11-04 21:34:46
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (40->22), mult. (35->25), div. (0->0), fcn. (48->10), ass. (0->17)
	t81 = sin(pkin(10));
	t86 = sin(qJ(1));
	t93 = t86 * t81;
	t83 = cos(pkin(10));
	t92 = t86 * t83;
	t88 = cos(qJ(1));
	t91 = t88 * t81;
	t90 = t88 * t83;
	t82 = sin(pkin(9));
	t84 = cos(pkin(9));
	t87 = cos(qJ(3));
	t89 = (pkin(4) * t84 + qJ(5) * t82 + pkin(3)) * sin(qJ(3)) - (-t82 * pkin(4) + qJ(5) * t84) * t87 + qJ(2);
	t80 = qJ(3) + pkin(9);
	t79 = pkin(1) + pkin(7) + qJ(4);
	t78 = cos(t80);
	t77 = sin(t80);
	t1 = [t77 * t92 + t91, -t77 * t93 + t90, -t86 * t78, t79 * t88 + t89 * t86 + 0; -t77 * t90 + t93, t77 * t91 + t92, t88 * t78, t79 * t86 - t89 * t88 + 0; t78 * t83, -t78 * t81, t77, t87 * pkin(3) + t78 * pkin(4) + t77 * qJ(5) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:34:46
	% EndTime: 2020-11-04 21:34:46
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (58->25), mult. (40->24), div. (0->0), fcn. (53->10), ass. (0->17)
	t105 = sin(qJ(1));
	t101 = pkin(10) + qJ(6);
	t96 = sin(t101);
	t112 = t105 * t96;
	t98 = cos(t101);
	t111 = t105 * t98;
	t106 = cos(qJ(1));
	t110 = t106 * t96;
	t109 = t106 * t98;
	t108 = sin(pkin(10)) * pkin(5) + pkin(1) + pkin(7) + qJ(4);
	t104 = -pkin(8) - qJ(5);
	t95 = cos(pkin(10)) * pkin(5) + pkin(4);
	t102 = qJ(3) + pkin(9);
	t97 = sin(t102);
	t99 = cos(t102);
	t107 = t104 * t99 + t95 * t97 + sin(qJ(3)) * pkin(3) + qJ(2);
	t1 = [t97 * t111 + t110, -t97 * t112 + t109, -t105 * t99, t107 * t105 + t108 * t106 + 0; -t97 * t109 + t112, t97 * t110 + t111, t106 * t99, t108 * t105 - t107 * t106 + 0; t99 * t98, -t99 * t96, t97, t99 * t95 - t97 * t104 + cos(qJ(3)) * pkin(3) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end