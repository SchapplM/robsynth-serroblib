% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPRP5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:01
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPPRP5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:01:51
	% EndTime: 2020-11-04 22:01:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:01:51
	% EndTime: 2020-11-04 22:01:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t54 = cos(qJ(1));
	t53 = sin(qJ(1));
	t1 = [t54, -t53, 0, 0; t53, t54, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:01:51
	% EndTime: 2020-11-04 22:01:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t58 = cos(qJ(1));
	t57 = cos(qJ(2));
	t56 = sin(qJ(1));
	t55 = sin(qJ(2));
	t1 = [t58 * t57, -t58 * t55, t56, t58 * pkin(1) + t56 * pkin(7) + 0; t56 * t57, -t56 * t55, -t58, t56 * pkin(1) - t58 * pkin(7) + 0; t55, t57, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:01:51
	% EndTime: 2020-11-04 22:01:51
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (16->14), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->6)
	t63 = cos(qJ(1));
	t62 = cos(qJ(2));
	t61 = sin(qJ(1));
	t60 = sin(qJ(2));
	t59 = t62 * pkin(2) + t60 * qJ(3) + pkin(1);
	t1 = [t61, -t63 * t62, t63 * t60, t61 * pkin(7) + t59 * t63 + 0; -t63, -t61 * t62, t61 * t60, -t63 * pkin(7) + t59 * t61 + 0; 0, -t60, -t62, t60 * pkin(2) - t62 * qJ(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:01:51
	% EndTime: 2020-11-04 22:01:51
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (22->17), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->14)
	t65 = sin(pkin(9));
	t69 = sin(qJ(1));
	t76 = t69 * t65;
	t66 = cos(pkin(9));
	t75 = t69 * t66;
	t71 = cos(qJ(1));
	t74 = t71 * t65;
	t73 = t71 * t66;
	t72 = pkin(3) + pkin(7);
	t70 = cos(qJ(2));
	t68 = sin(qJ(2));
	t67 = pkin(2) + qJ(4);
	t64 = t68 * qJ(3) + t67 * t70 + pkin(1);
	t1 = [t68 * t74 + t75, t68 * t73 - t76, t71 * t70, t64 * t71 + t72 * t69 + 0; t68 * t76 - t73, t68 * t75 + t74, t69 * t70, t64 * t69 - t72 * t71 + 0; -t70 * t65, -t70 * t66, t68, -t70 * qJ(3) + t67 * t68 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:01:51
	% EndTime: 2020-11-04 22:01:51
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (40->20), mult. (31->22), div. (0->0), fcn. (44->8), ass. (0->16)
	t83 = pkin(9) + qJ(5);
	t80 = sin(t83);
	t85 = sin(qJ(1));
	t91 = t85 * t80;
	t81 = cos(t83);
	t90 = t85 * t81;
	t87 = cos(qJ(1));
	t89 = t87 * t80;
	t88 = t87 * t81;
	t86 = cos(qJ(2));
	t84 = sin(qJ(2));
	t82 = pkin(2) + pkin(8) + qJ(4);
	t79 = sin(pkin(9)) * pkin(4) + qJ(3);
	t78 = cos(pkin(9)) * pkin(4) + pkin(3) + pkin(7);
	t77 = t79 * t84 + t82 * t86 + pkin(1);
	t1 = [t84 * t89 + t90, t84 * t88 - t91, t87 * t86, t77 * t87 + t78 * t85 + 0; t84 * t91 - t88, t84 * t90 + t89, t85 * t86, t77 * t85 - t78 * t87 + 0; -t86 * t80, -t86 * t81, t84, -t79 * t86 + t82 * t84 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:01:51
	% EndTime: 2020-11-04 22:01:51
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (59->26), mult. (51->28), div. (0->0), fcn. (68->8), ass. (0->20)
	t104 = sin(qJ(1));
	t102 = pkin(9) + qJ(5);
	t99 = sin(t102);
	t110 = t104 * t99;
	t106 = cos(qJ(1));
	t109 = t106 * t99;
	t100 = cos(t102);
	t108 = t104 * t100;
	t107 = t106 * t100;
	t105 = cos(qJ(2));
	t103 = sin(qJ(2));
	t101 = pkin(2) + pkin(8) + qJ(4);
	t98 = sin(pkin(9)) * pkin(4) + qJ(3);
	t97 = cos(pkin(9)) * pkin(4) + pkin(3) + pkin(7);
	t96 = t101 * t105 + t98 * t103 + pkin(1);
	t95 = t103 * t110 - t107;
	t94 = t103 * t108 + t109;
	t93 = t103 * t109 + t108;
	t92 = -t103 * t107 + t110;
	t1 = [t93, t106 * t105, t92, t93 * pkin(5) + t92 * qJ(6) + t97 * t104 + t96 * t106 + 0; t95, t104 * t105, -t94, t95 * pkin(5) - t94 * qJ(6) + t96 * t104 - t97 * t106 + 0; -t105 * t99, t103, t105 * t100, t101 * t103 + pkin(6) + 0 + (-pkin(5) * t99 + qJ(6) * t100 - t98) * t105; 0, 0, 0, 1;];
	Tc_mdh = t1;
end