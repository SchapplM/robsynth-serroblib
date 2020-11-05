% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPPR1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:22
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRPPR1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:22:22
	% EndTime: 2020-11-04 22:22:22
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:22:22
	% EndTime: 2020-11-04 22:22:22
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t59 = cos(qJ(1));
	t58 = sin(qJ(1));
	t1 = [t59, -t58, 0, 0; t58, t59, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:22:22
	% EndTime: 2020-11-04 22:22:22
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t63 = cos(qJ(1));
	t62 = cos(qJ(2));
	t61 = sin(qJ(1));
	t60 = sin(qJ(2));
	t1 = [t63 * t62, -t63 * t60, t61, t63 * pkin(1) + t61 * pkin(7) + 0; t61 * t62, -t61 * t60, -t63, t61 * pkin(1) - t63 * pkin(7) + 0; t60, t62, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:22:22
	% EndTime: 2020-11-04 22:22:22
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t70 = pkin(8) + pkin(7);
	t69 = cos(qJ(1));
	t68 = sin(qJ(1));
	t67 = qJ(2) + qJ(3);
	t66 = cos(t67);
	t65 = sin(t67);
	t64 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t69 * t66, -t69 * t65, t68, t69 * t64 + t70 * t68 + 0; t68 * t66, -t68 * t65, -t69, t68 * t64 - t69 * t70 + 0; t65, t66, 0, sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:22:22
	% EndTime: 2020-11-04 22:22:22
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->15), mult. (14->12), div. (0->0), fcn. (22->8), ass. (0->9)
	t76 = qJ(2) + qJ(3);
	t78 = cos(qJ(1));
	t77 = sin(qJ(1));
	t75 = -qJ(4) - pkin(8) - pkin(7);
	t74 = pkin(10) + t76;
	t73 = cos(t74);
	t72 = sin(t74);
	t71 = pkin(3) * cos(t76) + cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t78 * t73, -t78 * t72, t77, t78 * t71 - t77 * t75 + 0; t77 * t73, -t77 * t72, -t78, t77 * t71 + t78 * t75 + 0; t72, t73, 0, pkin(3) * sin(t76) + sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:22:22
	% EndTime: 2020-11-04 22:22:22
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (60->22), mult. (36->24), div. (0->0), fcn. (49->10), ass. (0->15)
	t85 = sin(pkin(11));
	t87 = sin(qJ(1));
	t93 = t87 * t85;
	t86 = cos(pkin(11));
	t92 = t87 * t86;
	t88 = cos(qJ(1));
	t91 = t88 * t85;
	t90 = t88 * t86;
	t84 = qJ(2) + qJ(3);
	t82 = pkin(10) + t84;
	t80 = sin(t82);
	t81 = cos(t82);
	t89 = pkin(4) * t81 + qJ(5) * t80 + pkin(3) * cos(t84) + cos(qJ(2)) * pkin(2) + pkin(1);
	t83 = -qJ(4) - pkin(8) - pkin(7);
	t1 = [t81 * t90 + t93, -t81 * t91 + t92, t88 * t80, -t87 * t83 + t89 * t88 + 0; t81 * t92 - t91, -t81 * t93 - t90, t87 * t80, t88 * t83 + t89 * t87 + 0; t80 * t86, -t80 * t85, -t81, t80 * pkin(4) - t81 * qJ(5) + pkin(3) * sin(t84) + sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:22:22
	% EndTime: 2020-11-04 22:22:22
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (78->27), mult. (43->26), div. (0->0), fcn. (56->12), ass. (0->18)
	t106 = sin(qJ(1));
	t102 = pkin(11) + qJ(6);
	t98 = sin(t102);
	t113 = t106 * t98;
	t99 = cos(t102);
	t112 = t106 * t99;
	t107 = cos(qJ(1));
	t111 = t107 * t98;
	t110 = t107 * t99;
	t103 = qJ(2) + qJ(3);
	t109 = sin(pkin(11)) * pkin(5) + pkin(7) + pkin(8) + qJ(4);
	t105 = -pkin(9) - qJ(5);
	t100 = pkin(10) + t103;
	t95 = sin(t100);
	t96 = cos(t100);
	t97 = cos(pkin(11)) * pkin(5) + pkin(4);
	t108 = -t105 * t95 + t96 * t97 + pkin(3) * cos(t103) + cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t96 * t110 + t113, -t96 * t111 + t112, t107 * t95, t109 * t106 + t108 * t107 + 0; t96 * t112 - t111, -t96 * t113 - t110, t106 * t95, t108 * t106 - t109 * t107 + 0; t95 * t99, -t95 * t98, -t96, t95 * t97 + t96 * t105 + pkin(3) * sin(t103) + sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end