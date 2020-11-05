% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPRP9 (for one body)
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

function Tc_mdh = S6RPRPRP9_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP9_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:38:29
	% EndTime: 2020-11-04 21:38:29
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:38:29
	% EndTime: 2020-11-04 21:38:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t52 = cos(qJ(1));
	t51 = sin(qJ(1));
	t1 = [t52, -t51, 0, 0; t51, t52, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:38:29
	% EndTime: 2020-11-04 21:38:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t54 = cos(qJ(1));
	t53 = sin(qJ(1));
	t1 = [0, -t54, t53, t54 * pkin(1) + t53 * qJ(2) + 0; 0, -t53, -t54, t53 * pkin(1) - t54 * qJ(2) + 0; 1, 0, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:38:29
	% EndTime: 2020-11-04 21:38:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->10), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->6)
	t59 = pkin(1) + pkin(7);
	t58 = cos(qJ(1));
	t57 = cos(qJ(3));
	t56 = sin(qJ(1));
	t55 = sin(qJ(3));
	t1 = [t56 * t55, t56 * t57, t58, t56 * qJ(2) + t59 * t58 + 0; -t58 * t55, -t58 * t57, t56, -t58 * qJ(2) + t59 * t56 + 0; t57, -t55, 0, pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:38:29
	% EndTime: 2020-11-04 21:38:29
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (20->17), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->13)
	t61 = sin(pkin(9));
	t64 = sin(qJ(1));
	t71 = t64 * t61;
	t62 = cos(pkin(9));
	t70 = t64 * t62;
	t66 = cos(qJ(1));
	t69 = t66 * t61;
	t68 = t66 * t62;
	t67 = pkin(1) + pkin(7);
	t65 = cos(qJ(3));
	t63 = sin(qJ(3));
	t60 = -t63 * pkin(3) + t65 * qJ(4) - qJ(2);
	t1 = [t63 * t70 + t69, -t63 * t71 + t68, -t64 * t65, -t60 * t64 + t67 * t66 + 0; -t63 * t68 + t71, t63 * t69 + t70, t66 * t65, t60 * t66 + t67 * t64 + 0; t65 * t62, -t65 * t61, t63, t65 * pkin(3) + t63 * qJ(4) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:38:29
	% EndTime: 2020-11-04 21:38:29
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (38->21), mult. (31->22), div. (0->0), fcn. (44->8), ass. (0->16)
	t76 = pkin(9) + qJ(5);
	t74 = sin(t76);
	t79 = sin(qJ(1));
	t86 = t79 * t74;
	t75 = cos(t76);
	t85 = t79 * t75;
	t81 = cos(qJ(1));
	t84 = t81 * t74;
	t83 = t81 * t75;
	t73 = cos(pkin(9)) * pkin(4) + pkin(3);
	t77 = pkin(8) + qJ(4);
	t78 = sin(qJ(3));
	t80 = cos(qJ(3));
	t82 = t73 * t78 - t77 * t80 + qJ(2);
	t72 = sin(pkin(9)) * pkin(4) + pkin(1) + pkin(7);
	t1 = [t78 * t85 + t84, -t78 * t86 + t83, -t79 * t80, t72 * t81 + t82 * t79 + 0; -t78 * t83 + t86, t78 * t84 + t85, t81 * t80, t72 * t79 - t82 * t81 + 0; t80 * t75, -t80 * t74, t78, t80 * t73 + t77 * t78 + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:38:29
	% EndTime: 2020-11-04 21:38:29
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (57->27), mult. (51->29), div. (0->0), fcn. (68->8), ass. (0->19)
	t95 = pkin(9) + qJ(5);
	t93 = sin(t95);
	t98 = sin(qJ(1));
	t104 = t98 * t93;
	t94 = cos(t95);
	t103 = t98 * t94;
	t100 = cos(qJ(1));
	t97 = sin(qJ(3));
	t102 = t100 * t97;
	t92 = cos(pkin(9)) * pkin(4) + pkin(3);
	t96 = pkin(8) + qJ(4);
	t99 = cos(qJ(3));
	t101 = t92 * t97 - t96 * t99 + qJ(2);
	t91 = sin(pkin(9)) * pkin(4) + pkin(1) + pkin(7);
	t90 = -t94 * t102 + t104;
	t89 = t93 * t102 + t103;
	t88 = t100 * t93 + t97 * t103;
	t87 = -t100 * t94 + t97 * t104;
	t1 = [t88, -t98 * t99, t87, t88 * pkin(5) + t87 * qJ(6) + t91 * t100 + t101 * t98 + 0; t90, t100 * t99, -t89, t90 * pkin(5) - t89 * qJ(6) - t101 * t100 + t91 * t98 + 0; t99 * t94, t97, t99 * t93, t96 * t97 + pkin(2) + pkin(6) + 0 + (pkin(5) * t94 + qJ(6) * t93 + t92) * t99; 0, 0, 0, 1;];
	Tc_mdh = t1;
end