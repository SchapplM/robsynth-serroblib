% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPRP8 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:38
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRPRP8_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP8_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:38:09
	% EndTime: 2020-11-04 21:38:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:38:09
	% EndTime: 2020-11-04 21:38:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t58 = cos(qJ(1));
	t57 = sin(qJ(1));
	t1 = [t58, -t57, 0, 0; t57, t58, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:38:09
	% EndTime: 2020-11-04 21:38:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t60 = cos(qJ(1));
	t59 = sin(qJ(1));
	t1 = [0, -t60, t59, t60 * pkin(1) + t59 * qJ(2) + 0; 0, -t59, -t60, t59 * pkin(1) - t60 * qJ(2) + 0; 1, 0, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:38:09
	% EndTime: 2020-11-04 21:38:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->10), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->6)
	t65 = pkin(1) + pkin(7);
	t64 = cos(qJ(1));
	t63 = cos(qJ(3));
	t62 = sin(qJ(1));
	t61 = sin(qJ(3));
	t1 = [t62 * t61, t62 * t63, t64, t62 * qJ(2) + t65 * t64 + 0; -t64 * t61, -t64 * t63, t62, -t64 * qJ(2) + t65 * t62 + 0; t63, -t61, 0, pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:38:09
	% EndTime: 2020-11-04 21:38:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->13), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t72 = cos(qJ(1));
	t71 = sin(qJ(1));
	t70 = qJ(3) + pkin(9);
	t69 = pkin(1) + pkin(7) + qJ(4);
	t68 = cos(t70);
	t67 = sin(t70);
	t66 = sin(qJ(3)) * pkin(3) + qJ(2);
	t1 = [t71 * t67, t71 * t68, t72, t66 * t71 + t69 * t72 + 0; -t72 * t67, -t72 * t68, t71, -t66 * t72 + t69 * t71 + 0; t68, -t67, 0, cos(qJ(3)) * pkin(3) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:38:09
	% EndTime: 2020-11-04 21:38:09
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (40->22), mult. (35->25), div. (0->0), fcn. (48->10), ass. (0->17)
	t81 = sin(qJ(5));
	t83 = sin(qJ(1));
	t91 = t83 * t81;
	t84 = cos(qJ(5));
	t90 = t83 * t84;
	t86 = cos(qJ(1));
	t89 = t86 * t81;
	t88 = t86 * t84;
	t79 = sin(pkin(9));
	t80 = cos(pkin(9));
	t85 = cos(qJ(3));
	t87 = (t80 * pkin(4) + t79 * pkin(8) + pkin(3)) * sin(qJ(3)) - (-t79 * pkin(4) + t80 * pkin(8)) * t85 + qJ(2);
	t78 = qJ(3) + pkin(9);
	t77 = pkin(1) + pkin(7) + qJ(4);
	t76 = cos(t78);
	t75 = sin(t78);
	t1 = [t75 * t90 + t89, -t75 * t91 + t88, -t83 * t76, t77 * t86 + t87 * t83 + 0; -t75 * t88 + t91, t75 * t89 + t90, t86 * t76, t77 * t83 - t87 * t86 + 0; t76 * t84, -t76 * t81, t75, t85 * pkin(3) + t76 * pkin(4) + t75 * pkin(8) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:38:09
	% EndTime: 2020-11-04 21:38:09
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (55->28), mult. (55->31), div. (0->0), fcn. (72->10), ass. (0->21)
	t104 = sin(qJ(5));
	t106 = sin(qJ(1));
	t114 = t106 * t104;
	t107 = cos(qJ(5));
	t113 = t106 * t107;
	t109 = cos(qJ(1));
	t112 = t109 * t104;
	t111 = t109 * t107;
	t102 = sin(pkin(9));
	t103 = cos(pkin(9));
	t108 = cos(qJ(3));
	t110 = (t103 * pkin(4) + t102 * pkin(8) + pkin(3)) * sin(qJ(3)) - (-t102 * pkin(4) + t103 * pkin(8)) * t108 + qJ(2);
	t101 = qJ(3) + pkin(9);
	t100 = pkin(1) + pkin(7) + qJ(4);
	t99 = cos(t101);
	t98 = sin(t101);
	t95 = -t98 * t111 + t114;
	t94 = t98 * t112 + t113;
	t93 = t98 * t113 + t112;
	t92 = t98 * t114 - t111;
	t1 = [t93, -t106 * t99, t92, t93 * pkin(5) + t92 * qJ(6) + t100 * t109 + t110 * t106 + 0; t95, t109 * t99, -t94, t95 * pkin(5) - t94 * qJ(6) + t100 * t106 - t110 * t109 + 0; t99 * t107, t98, t99 * t104, t108 * pkin(3) + t98 * pkin(8) + pkin(2) + pkin(6) + 0 + (pkin(5) * t107 + qJ(6) * t104 + pkin(4)) * t99; 0, 0, 0, 1;];
	Tc_mdh = t1;
end