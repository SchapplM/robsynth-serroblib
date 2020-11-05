% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PPRPR3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:53
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PPRPR3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:53:05
	% EndTime: 2020-11-04 19:53:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:53:05
	% EndTime: 2020-11-04 19:53:05
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t56 = cos(pkin(7));
	t55 = sin(pkin(7));
	t1 = [t56, -t55, 0, 0; t55, t56, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:53:05
	% EndTime: 2020-11-04 19:53:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t60 = cos(pkin(7));
	t59 = cos(pkin(8));
	t58 = sin(pkin(7));
	t57 = sin(pkin(8));
	t1 = [t60 * t59, -t60 * t57, t58, t60 * pkin(1) + t58 * qJ(2) + 0; t58 * t59, -t58 * t57, -t60, t58 * pkin(1) - t60 * qJ(2) + 0; t57, t59, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:53:05
	% EndTime: 2020-11-04 19:53:05
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->15), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->12)
	t63 = sin(pkin(7));
	t66 = sin(qJ(3));
	t71 = t63 * t66;
	t67 = cos(qJ(3));
	t70 = t63 * t67;
	t65 = cos(pkin(7));
	t69 = t65 * t66;
	t68 = t65 * t67;
	t64 = cos(pkin(8));
	t62 = sin(pkin(8));
	t61 = t64 * pkin(2) + t62 * pkin(5) + pkin(1);
	t1 = [t64 * t68 + t71, -t64 * t69 + t70, t65 * t62, t63 * qJ(2) + t61 * t65 + 0; t64 * t70 - t69, -t64 * t71 - t68, t63 * t62, -t65 * qJ(2) + t61 * t63 + 0; t62 * t67, -t62 * t66, -t64, t62 * pkin(2) - t64 * pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:53:05
	% EndTime: 2020-11-04 19:53:05
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (35->19), mult. (37->23), div. (0->0), fcn. (50->8), ass. (0->15)
	t77 = sin(pkin(7));
	t78 = cos(pkin(8));
	t86 = t77 * t78;
	t75 = qJ(3) + pkin(9);
	t73 = sin(t75);
	t79 = cos(pkin(7));
	t85 = t79 * t73;
	t74 = cos(t75);
	t84 = t79 * t74;
	t83 = pkin(3) * sin(qJ(3)) + qJ(2);
	t72 = cos(qJ(3)) * pkin(3) + pkin(2);
	t76 = sin(pkin(8));
	t80 = -qJ(4) - pkin(5);
	t82 = t72 * t78 - t76 * t80 + pkin(1);
	t1 = [t77 * t73 + t78 * t84, t77 * t74 - t78 * t85, t79 * t76, t83 * t77 + t82 * t79 + 0; t74 * t86 - t85, -t73 * t86 - t84, t77 * t76, t82 * t77 - t83 * t79 + 0; t76 * t74, -t76 * t73, -t78, t76 * t72 + t78 * t80 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:53:05
	% EndTime: 2020-11-04 19:53:05
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (66->29), mult. (78->40), div. (0->0), fcn. (103->10), ass. (0->23)
	t96 = sin(pkin(7));
	t97 = cos(pkin(8));
	t109 = t96 * t97;
	t94 = qJ(3) + pkin(9);
	t92 = sin(t94);
	t98 = cos(pkin(7));
	t108 = t98 * t92;
	t93 = cos(t94);
	t107 = t98 * t93;
	t100 = sin(qJ(5));
	t95 = sin(pkin(8));
	t106 = t100 * t95;
	t102 = cos(qJ(5));
	t105 = t102 * t95;
	t104 = pkin(3) * sin(qJ(3)) + qJ(2);
	t91 = cos(qJ(3)) * pkin(3) + pkin(2);
	t99 = -qJ(4) - pkin(5);
	t103 = t91 * t97 - t95 * t99 + pkin(1);
	t90 = t97 * t107 + t96 * t92;
	t89 = t97 * t108 - t96 * t93;
	t88 = t93 * t109 - t108;
	t87 = t92 * t109 + t107;
	t1 = [t90 * t102 + t98 * t106, -t90 * t100 + t98 * t105, t89, t90 * pkin(4) + t89 * pkin(6) + t103 * t98 + t104 * t96 + 0; t88 * t102 + t96 * t106, -t88 * t100 + t96 * t105, t87, t88 * pkin(4) + t87 * pkin(6) + t103 * t96 - t104 * t98 + 0; -t97 * t100 + t93 * t105, -t97 * t102 - t93 * t106, t95 * t92, t97 * t99 + qJ(1) + 0 + (pkin(4) * t93 + pkin(6) * t92 + t91) * t95; 0, 0, 0, 1;];
	Tc_mdh = t1;
end