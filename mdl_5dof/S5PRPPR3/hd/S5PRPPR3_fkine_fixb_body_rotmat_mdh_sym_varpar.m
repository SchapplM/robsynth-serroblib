% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRPPR3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:56
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PRPPR3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:56:58
	% EndTime: 2020-11-04 19:56:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:56:58
	% EndTime: 2020-11-04 19:56:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t50 = cos(pkin(7));
	t49 = sin(pkin(7));
	t1 = [t50, -t49, 0, 0; t49, t50, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:56:58
	% EndTime: 2020-11-04 19:56:58
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t54 = cos(qJ(2));
	t53 = sin(qJ(2));
	t52 = cos(pkin(7));
	t51 = sin(pkin(7));
	t1 = [t52 * t54, -t52 * t53, t51, t52 * pkin(1) + t51 * pkin(5) + 0; t51 * t54, -t51 * t53, -t52, t51 * pkin(1) - t52 * pkin(5) + 0; t53, t54, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:56:58
	% EndTime: 2020-11-04 19:56:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t61 = -qJ(3) - pkin(5);
	t60 = cos(pkin(7));
	t59 = sin(pkin(7));
	t58 = qJ(2) + pkin(8);
	t57 = cos(t58);
	t56 = sin(t58);
	t55 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t60 * t57, -t60 * t56, t59, t60 * t55 - t59 * t61 + 0; t59 * t57, -t59 * t56, -t60, t59 * t55 + t60 * t61 + 0; t56, t57, 0, sin(qJ(2)) * pkin(2) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:56:58
	% EndTime: 2020-11-04 19:56:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->18), mult. (21->14), div. (0->0), fcn. (29->6), ass. (0->8)
	t65 = qJ(2) + pkin(8);
	t63 = sin(t65);
	t64 = cos(t65);
	t69 = pkin(3) * t64 + qJ(4) * t63 + cos(qJ(2)) * pkin(2) + pkin(1);
	t68 = -qJ(3) - pkin(5);
	t67 = cos(pkin(7));
	t66 = sin(pkin(7));
	t1 = [t66, -t67 * t64, t67 * t63, -t66 * t68 + t69 * t67 + 0; -t67, -t66 * t64, t66 * t63, t69 * t66 + t67 * t68 + 0; 0, -t63, -t64, t63 * pkin(3) - t64 * qJ(4) + sin(qJ(2)) * pkin(2) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:56:58
	% EndTime: 2020-11-04 19:56:58
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (44->22), mult. (39->25), div. (0->0), fcn. (52->10), ass. (0->18)
	t77 = sin(pkin(7));
	t80 = sin(qJ(5));
	t89 = t77 * t80;
	t82 = cos(qJ(5));
	t88 = t77 * t82;
	t79 = cos(pkin(7));
	t87 = t79 * t80;
	t86 = t79 * t82;
	t76 = sin(pkin(8));
	t78 = cos(pkin(8));
	t81 = sin(qJ(2));
	t84 = pkin(3) + pkin(6);
	t85 = (qJ(4) * t76 + t84 * t78 + pkin(2)) * cos(qJ(2)) + (qJ(4) * t78 - t76 * t84) * t81 + pkin(1);
	t75 = qJ(2) + pkin(8);
	t74 = qJ(3) + pkin(4) + pkin(5);
	t73 = cos(t75);
	t72 = sin(t75);
	t1 = [t72 * t87 + t88, t72 * t86 - t89, t79 * t73, t74 * t77 + t85 * t79 + 0; t72 * t89 - t86, t72 * t88 + t87, t77 * t73, -t74 * t79 + t85 * t77 + 0; -t73 * t80, -t73 * t82, t72, t81 * pkin(2) - t73 * qJ(4) + t84 * t72 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end