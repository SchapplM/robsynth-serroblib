% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRPRP6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:59
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PRPRP6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:59:09
	% EndTime: 2020-11-04 19:59:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:59:09
	% EndTime: 2020-11-04 19:59:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t51 = cos(pkin(7));
	t50 = sin(pkin(7));
	t1 = [t51, -t50, 0, 0; t50, t51, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:59:09
	% EndTime: 2020-11-04 19:59:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t55 = cos(qJ(2));
	t54 = sin(qJ(2));
	t53 = cos(pkin(7));
	t52 = sin(pkin(7));
	t1 = [t53 * t55, -t53 * t54, t52, t53 * pkin(1) + t52 * pkin(5) + 0; t52 * t55, -t52 * t54, -t53, t52 * pkin(1) - t53 * pkin(5) + 0; t54, t55, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:59:09
	% EndTime: 2020-11-04 19:59:09
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (16->14), mult. (18->12), div. (0->0), fcn. (26->4), ass. (0->6)
	t58 = sin(qJ(2));
	t59 = cos(qJ(2));
	t60 = pkin(2) * t59 + qJ(3) * t58 + pkin(1);
	t57 = cos(pkin(7));
	t56 = sin(pkin(7));
	t1 = [t56, -t57 * t59, t57 * t58, t56 * pkin(5) + t60 * t57 + 0; -t57, -t56 * t59, t56 * t58, -t57 * pkin(5) + t60 * t56 + 0; 0, -t58, -t59, t58 * pkin(2) - t59 * qJ(3) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:59:09
	% EndTime: 2020-11-04 19:59:09
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (22->17), mult. (30->22), div. (0->0), fcn. (43->6), ass. (0->12)
	t63 = sin(qJ(4));
	t64 = sin(qJ(2));
	t71 = t63 * t64;
	t65 = cos(qJ(4));
	t70 = t64 * t65;
	t66 = cos(qJ(2));
	t68 = pkin(2) + pkin(6);
	t69 = qJ(3) * t64 + t66 * t68 + pkin(1);
	t67 = pkin(3) + pkin(5);
	t62 = cos(pkin(7));
	t61 = sin(pkin(7));
	t1 = [t61 * t65 + t62 * t71, -t61 * t63 + t62 * t70, t62 * t66, t67 * t61 + t69 * t62 + 0; t61 * t71 - t62 * t65, t61 * t70 + t62 * t63, t61 * t66, t69 * t61 - t67 * t62 + 0; -t66 * t63, -t66 * t65, t64, -t66 * qJ(3) + t68 * t64 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:59:09
	% EndTime: 2020-11-04 19:59:09
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (31->20), mult. (44->26), div. (0->0), fcn. (57->6), ass. (0->13)
	t76 = sin(qJ(4));
	t77 = sin(qJ(2));
	t85 = t76 * t77;
	t78 = cos(qJ(4));
	t84 = t77 * t78;
	t72 = -t76 * pkin(4) + t78 * qJ(5) - qJ(3);
	t79 = cos(qJ(2));
	t81 = pkin(2) + pkin(6);
	t83 = -t72 * t77 + t79 * t81 + pkin(1);
	t82 = pkin(4) * t78 + qJ(5) * t76 + pkin(3) + pkin(5);
	t75 = cos(pkin(7));
	t74 = sin(pkin(7));
	t1 = [t74 * t78 + t75 * t85, t75 * t79, t74 * t76 - t75 * t84, t82 * t74 + t83 * t75 + 0; t74 * t85 - t75 * t78, t74 * t79, -t74 * t84 - t75 * t76, t83 * t74 - t82 * t75 + 0; -t79 * t76, t77, t79 * t78, t72 * t79 + t81 * t77 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end