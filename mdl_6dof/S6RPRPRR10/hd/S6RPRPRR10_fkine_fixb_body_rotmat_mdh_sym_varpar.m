% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPRR10 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:42
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRPRR10_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR10_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:42:20
	% EndTime: 2020-11-04 21:42:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:42:20
	% EndTime: 2020-11-04 21:42:20
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t48 = cos(qJ(1));
	t47 = sin(qJ(1));
	t1 = [t48, -t47, 0, 0; t47, t48, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:42:20
	% EndTime: 2020-11-04 21:42:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t50 = cos(qJ(1));
	t49 = sin(qJ(1));
	t1 = [0, -t50, t49, t50 * pkin(1) + t49 * qJ(2) + 0; 0, -t49, -t50, t49 * pkin(1) - t50 * qJ(2) + 0; 1, 0, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:42:20
	% EndTime: 2020-11-04 21:42:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->10), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->6)
	t55 = pkin(1) + pkin(7);
	t54 = cos(qJ(1));
	t53 = cos(qJ(3));
	t52 = sin(qJ(1));
	t51 = sin(qJ(3));
	t1 = [t52 * t51, t52 * t53, t54, t52 * qJ(2) + t55 * t54 + 0; -t54 * t51, -t54 * t53, t52, -t54 * qJ(2) + t55 * t52 + 0; t53, -t51, 0, pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:42:20
	% EndTime: 2020-11-04 21:42:20
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (20->17), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->13)
	t57 = sin(pkin(10));
	t60 = sin(qJ(1));
	t67 = t60 * t57;
	t58 = cos(pkin(10));
	t66 = t60 * t58;
	t62 = cos(qJ(1));
	t65 = t62 * t57;
	t64 = t62 * t58;
	t63 = pkin(1) + pkin(7);
	t61 = cos(qJ(3));
	t59 = sin(qJ(3));
	t56 = -t59 * pkin(3) + t61 * qJ(4) - qJ(2);
	t1 = [t59 * t66 + t65, -t59 * t67 + t64, -t60 * t61, -t56 * t60 + t63 * t62 + 0; -t59 * t64 + t67, t59 * t65 + t66, t62 * t61, t56 * t62 + t63 * t60 + 0; t61 * t58, -t61 * t57, t59, t61 * pkin(3) + t59 * qJ(4) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:42:20
	% EndTime: 2020-11-04 21:42:20
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (38->21), mult. (31->22), div. (0->0), fcn. (44->8), ass. (0->16)
	t72 = pkin(10) + qJ(5);
	t70 = sin(t72);
	t75 = sin(qJ(1));
	t82 = t75 * t70;
	t71 = cos(t72);
	t81 = t75 * t71;
	t77 = cos(qJ(1));
	t80 = t77 * t70;
	t79 = t77 * t71;
	t69 = cos(pkin(10)) * pkin(4) + pkin(3);
	t73 = pkin(8) + qJ(4);
	t74 = sin(qJ(3));
	t76 = cos(qJ(3));
	t78 = t69 * t74 - t73 * t76 + qJ(2);
	t68 = sin(pkin(10)) * pkin(4) + pkin(1) + pkin(7);
	t1 = [t74 * t81 + t80, -t74 * t82 + t79, -t75 * t76, t68 * t77 + t78 * t75 + 0; -t74 * t79 + t82, t74 * t80 + t81, t77 * t76, t68 * t75 - t78 * t77 + 0; t76 * t71, -t76 * t70, t74, t76 * t69 + t73 * t74 + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:42:20
	% EndTime: 2020-11-04 21:42:20
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (61->24), mult. (42->24), div. (0->0), fcn. (55->10), ass. (0->17)
	t89 = pkin(10) + qJ(5);
	t87 = qJ(6) + t89;
	t85 = sin(t87);
	t91 = sin(qJ(1));
	t100 = t91 * t85;
	t86 = cos(t87);
	t99 = t91 * t86;
	t93 = cos(qJ(1));
	t98 = t93 * t85;
	t97 = t93 * t86;
	t96 = pkin(5) * sin(t89) + sin(pkin(10)) * pkin(4) + pkin(1) + pkin(7);
	t83 = pkin(5) * cos(t89) + cos(pkin(10)) * pkin(4) + pkin(3);
	t88 = -pkin(9) - pkin(8) - qJ(4);
	t90 = sin(qJ(3));
	t92 = cos(qJ(3));
	t95 = t83 * t90 + t88 * t92 + qJ(2);
	t1 = [t90 * t99 + t98, -t90 * t100 + t97, -t91 * t92, t95 * t91 + t96 * t93 + 0; -t90 * t97 + t100, t90 * t98 + t99, t93 * t92, t96 * t91 - t95 * t93 + 0; t92 * t86, -t92 * t85, t90, t92 * t83 - t90 * t88 + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end