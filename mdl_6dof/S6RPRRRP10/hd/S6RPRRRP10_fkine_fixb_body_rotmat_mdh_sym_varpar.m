% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRRP10 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:53
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRRRP10_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP10_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:53:48
	% EndTime: 2020-11-04 21:53:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:53:48
	% EndTime: 2020-11-04 21:53:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t54 = cos(qJ(1));
	t53 = sin(qJ(1));
	t1 = [t54, -t53, 0, 0; t53, t54, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:53:48
	% EndTime: 2020-11-04 21:53:48
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t56 = cos(qJ(1));
	t55 = sin(qJ(1));
	t1 = [0, -t56, t55, pkin(1) * t56 + qJ(2) * t55 + 0; 0, -t55, -t56, pkin(1) * t55 - qJ(2) * t56 + 0; 1, 0, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:53:48
	% EndTime: 2020-11-04 21:53:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->10), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->6)
	t61 = pkin(1) + pkin(7);
	t60 = cos(qJ(1));
	t59 = cos(qJ(3));
	t58 = sin(qJ(1));
	t57 = sin(qJ(3));
	t1 = [t58 * t57, t58 * t59, t60, t58 * qJ(2) + t61 * t60 + 0; -t60 * t57, -t60 * t59, t58, -t60 * qJ(2) + t61 * t58 + 0; t59, -t57, 0, pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:53:48
	% EndTime: 2020-11-04 21:53:48
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (20->17), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->13)
	t63 = sin(qJ(4));
	t65 = sin(qJ(1));
	t73 = t65 * t63;
	t66 = cos(qJ(4));
	t72 = t65 * t66;
	t68 = cos(qJ(1));
	t71 = t68 * t63;
	t70 = t68 * t66;
	t69 = pkin(1) + pkin(7);
	t67 = cos(qJ(3));
	t64 = sin(qJ(3));
	t62 = -t64 * pkin(3) + t67 * pkin(8) - qJ(2);
	t1 = [t64 * t72 + t71, -t64 * t73 + t70, -t65 * t67, -t62 * t65 + t69 * t68 + 0; -t64 * t70 + t73, t64 * t71 + t72, t68 * t67, t62 * t68 + t69 * t65 + 0; t67 * t66, -t67 * t63, t64, t67 * pkin(3) + t64 * pkin(8) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:53:48
	% EndTime: 2020-11-04 21:53:48
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (38->21), mult. (31->22), div. (0->0), fcn. (44->8), ass. (0->16)
	t78 = qJ(4) + qJ(5);
	t76 = sin(t78);
	t80 = sin(qJ(1));
	t88 = t80 * t76;
	t77 = cos(t78);
	t87 = t80 * t77;
	t82 = cos(qJ(1));
	t86 = t82 * t76;
	t85 = t82 * t77;
	t75 = cos(qJ(4)) * pkin(4) + pkin(3);
	t79 = sin(qJ(3));
	t81 = cos(qJ(3));
	t83 = pkin(9) + pkin(8);
	t84 = t75 * t79 - t83 * t81 + qJ(2);
	t74 = sin(qJ(4)) * pkin(4) + pkin(1) + pkin(7);
	t1 = [t79 * t87 + t86, -t79 * t88 + t85, -t80 * t81, t74 * t82 + t84 * t80 + 0; -t79 * t85 + t88, t79 * t86 + t87, t82 * t81, t74 * t80 - t84 * t82 + 0; t81 * t77, -t81 * t76, t79, t81 * t75 + t83 * t79 + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:53:48
	% EndTime: 2020-11-04 21:53:48
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (57->27), mult. (51->29), div. (0->0), fcn. (68->8), ass. (0->19)
	t97 = qJ(4) + qJ(5);
	t95 = sin(t97);
	t99 = sin(qJ(1));
	t106 = t99 * t95;
	t96 = cos(t97);
	t105 = t99 * t96;
	t101 = cos(qJ(1));
	t98 = sin(qJ(3));
	t104 = t101 * t98;
	t100 = cos(qJ(3));
	t102 = pkin(9) + pkin(8);
	t94 = cos(qJ(4)) * pkin(4) + pkin(3);
	t103 = t102 * t100 - t94 * t98 - qJ(2);
	t93 = sin(qJ(4)) * pkin(4) + pkin(1) + pkin(7);
	t92 = -t96 * t104 + t106;
	t91 = t95 * t104 + t105;
	t90 = t101 * t95 + t98 * t105;
	t89 = -t101 * t96 + t98 * t106;
	t1 = [t90, -t99 * t100, t89, t90 * pkin(5) + t89 * qJ(6) + t93 * t101 - t103 * t99 + 0; t92, t101 * t100, -t91, t92 * pkin(5) - t91 * qJ(6) + t103 * t101 + t93 * t99 + 0; t100 * t96, t98, t100 * t95, t102 * t98 + pkin(2) + pkin(6) + 0 + (pkin(5) * t96 + qJ(6) * t95 + t94) * t100; 0, 0, 0, 1;];
	Tc_mdh = t1;
end