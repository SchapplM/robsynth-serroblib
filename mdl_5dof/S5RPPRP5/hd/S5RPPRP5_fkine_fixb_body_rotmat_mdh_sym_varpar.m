% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPRP5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:13
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPPRP5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:13:09
	% EndTime: 2020-11-04 20:13:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:13:09
	% EndTime: 2020-11-04 20:13:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t48 = cos(qJ(1));
	t47 = sin(qJ(1));
	t1 = [t48, -t47, 0, 0; t47, t48, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:13:09
	% EndTime: 2020-11-04 20:13:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t52 = cos(qJ(1));
	t51 = sin(qJ(1));
	t50 = cos(pkin(7));
	t49 = sin(pkin(7));
	t1 = [t52 * t50, -t52 * t49, t51, t52 * pkin(1) + t51 * qJ(2) + 0; t51 * t50, -t51 * t49, -t52, t51 * pkin(1) - t52 * qJ(2) + 0; t49, t50, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:13:09
	% EndTime: 2020-11-04 20:13:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (13->11), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->6)
	t57 = cos(qJ(1));
	t56 = sin(qJ(1));
	t55 = cos(pkin(7));
	t54 = sin(pkin(7));
	t53 = pkin(2) * t55 + t54 * qJ(3) + pkin(1);
	t1 = [t57 * t55, t56, t57 * t54, t56 * qJ(2) + t53 * t57 + 0; t56 * t55, -t57, t56 * t54, -t57 * qJ(2) + t53 * t56 + 0; t54, 0, -t55, t54 * pkin(2) - t55 * qJ(3) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:13:09
	% EndTime: 2020-11-04 20:13:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (23->16), mult. (26->16), div. (0->0), fcn. (40->6), ass. (0->12)
	t70 = cos(qJ(4));
	t69 = sin(qJ(4));
	t68 = pkin(2) + pkin(3);
	t67 = cos(qJ(1));
	t66 = sin(qJ(1));
	t65 = pkin(6) - qJ(2);
	t64 = cos(pkin(7));
	t63 = sin(pkin(7));
	t60 = t63 * qJ(3) + t68 * t64 + pkin(1);
	t59 = t63 * t70 - t64 * t69;
	t58 = -t63 * t69 - t64 * t70;
	t1 = [-t67 * t58, t67 * t59, -t66, t60 * t67 - t65 * t66 + 0; -t66 * t58, t66 * t59, t67, t60 * t66 + t65 * t67 + 0; t59, t58, 0, -t64 * qJ(3) + t68 * t63 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:13:09
	% EndTime: 2020-11-04 20:13:09
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (35->22), mult. (44->24), div. (0->0), fcn. (58->6), ass. (0->14)
	t84 = pkin(2) + pkin(3);
	t83 = cos(qJ(1));
	t82 = cos(qJ(4));
	t81 = sin(qJ(1));
	t80 = sin(qJ(4));
	t79 = pkin(6) - qJ(2);
	t78 = cos(pkin(7));
	t77 = sin(pkin(7));
	t75 = pkin(4) * t77 + t78 * qJ(5);
	t74 = t78 * pkin(4) - t77 * qJ(5);
	t73 = t77 * t80 + t78 * t82;
	t72 = t77 * t82 - t78 * t80;
	t71 = t77 * qJ(3) + t74 * t82 + t75 * t80 + t84 * t78 + pkin(1);
	t1 = [t83 * t73, -t81, -t83 * t72, t71 * t83 - t79 * t81 + 0; t81 * t73, t83, -t81 * t72, t71 * t81 + t79 * t83 + 0; t72, 0, t73, -t78 * qJ(3) - t74 * t80 + t75 * t82 + t84 * t77 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end