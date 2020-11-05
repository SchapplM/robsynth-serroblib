% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPPR3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:23
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRPPR3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:23:03
	% EndTime: 2020-11-04 22:23:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:23:03
	% EndTime: 2020-11-04 22:23:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t52 = cos(qJ(1));
	t51 = sin(qJ(1));
	t1 = [t52, -t51, 0, 0; t51, t52, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:23:03
	% EndTime: 2020-11-04 22:23:03
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t56 = cos(qJ(1));
	t55 = cos(qJ(2));
	t54 = sin(qJ(1));
	t53 = sin(qJ(2));
	t1 = [t56 * t55, -t56 * t53, t54, t56 * pkin(1) + t54 * pkin(7) + 0; t54 * t55, -t54 * t53, -t56, t54 * pkin(1) - t56 * pkin(7) + 0; t53, t55, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:23:03
	% EndTime: 2020-11-04 22:23:03
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t63 = pkin(8) + pkin(7);
	t62 = cos(qJ(1));
	t61 = sin(qJ(1));
	t60 = qJ(2) + qJ(3);
	t59 = cos(t60);
	t58 = sin(t60);
	t57 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t62 * t59, -t62 * t58, t61, t62 * t57 + t63 * t61 + 0; t61 * t59, -t61 * t58, -t62, t61 * t57 - t62 * t63 + 0; t58, t59, 0, sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:23:03
	% EndTime: 2020-11-04 22:23:03
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (30->15), mult. (21->14), div. (0->0), fcn. (29->6), ass. (0->8)
	t67 = qJ(2) + qJ(3);
	t65 = sin(t67);
	t66 = cos(t67);
	t71 = pkin(3) * t66 + qJ(4) * t65 + cos(qJ(2)) * pkin(2) + pkin(1);
	t70 = pkin(8) + pkin(7);
	t69 = cos(qJ(1));
	t68 = sin(qJ(1));
	t1 = [t69 * t66, t68, t69 * t65, t70 * t68 + t71 * t69 + 0; t68 * t66, -t69, t68 * t65, t71 * t68 - t69 * t70 + 0; t65, 0, -t66, t65 * pkin(3) - t66 * qJ(4) + sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:23:03
	% EndTime: 2020-11-04 22:23:03
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (40->21), mult. (25->17), div. (0->0), fcn. (33->8), ass. (0->12)
	t82 = pkin(3) + pkin(4);
	t81 = cos(qJ(1));
	t80 = cos(qJ(3));
	t79 = sin(qJ(1));
	t78 = sin(qJ(2));
	t77 = sin(qJ(3));
	t76 = qJ(2) + qJ(3);
	t75 = qJ(5) - pkin(7) - pkin(8);
	t74 = cos(t76);
	t73 = sin(t76);
	t72 = (qJ(4) * t77 + t82 * t80 + pkin(2)) * cos(qJ(2)) + pkin(1) + (qJ(4) * t80 - t77 * t82) * t78;
	t1 = [t81 * t73, -t81 * t74, -t79, t72 * t81 - t75 * t79 + 0; t79 * t73, -t79 * t74, t81, t72 * t79 + t75 * t81 + 0; -t74, -t73, 0, t78 * pkin(2) - t74 * qJ(4) + t82 * t73 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:23:03
	% EndTime: 2020-11-04 22:23:03
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (53->22), mult. (37->25), div. (0->0), fcn. (50->10), ass. (0->19)
	t90 = sin(qJ(6));
	t93 = sin(qJ(1));
	t100 = t93 * t90;
	t94 = cos(qJ(6));
	t99 = t93 * t94;
	t96 = cos(qJ(1));
	t98 = t96 * t90;
	t97 = t96 * t94;
	t95 = cos(qJ(3));
	t92 = sin(qJ(2));
	t91 = sin(qJ(3));
	t89 = qJ(4) + pkin(5);
	t88 = qJ(2) + qJ(3);
	t87 = pkin(3) + pkin(4) + pkin(9);
	t86 = qJ(5) - pkin(7) - pkin(8);
	t85 = cos(t88);
	t84 = sin(t88);
	t83 = (t87 * t95 + t89 * t91 + pkin(2)) * cos(qJ(2)) + pkin(1) + (-t91 * t87 + t89 * t95) * t92;
	t1 = [t84 * t97 - t100, -t84 * t98 - t99, t96 * t85, t83 * t96 - t86 * t93 + 0; t84 * t99 + t98, -t84 * t100 + t97, t93 * t85, t83 * t93 + t86 * t96 + 0; -t85 * t94, t85 * t90, t84, t92 * pkin(2) + t87 * t84 - t89 * t85 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end