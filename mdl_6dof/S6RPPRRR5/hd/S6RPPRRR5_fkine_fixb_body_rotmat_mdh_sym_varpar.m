% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPRRR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:31
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPPRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:31:41
	% EndTime: 2020-11-04 21:31:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:31:41
	% EndTime: 2020-11-04 21:31:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t46 = cos(qJ(1));
	t45 = sin(qJ(1));
	t1 = [t46, -t45, 0, 0; t45, t46, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:31:41
	% EndTime: 2020-11-04 21:31:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t48 = cos(qJ(1));
	t47 = sin(qJ(1));
	t1 = [0, -t48, t47, t48 * pkin(1) + t47 * qJ(2) + 0; 0, -t47, -t48, t47 * pkin(1) - t48 * qJ(2) + 0; 1, 0, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:31:41
	% EndTime: 2020-11-04 21:31:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->4)
	t51 = cos(qJ(1));
	t50 = sin(qJ(1));
	t49 = pkin(1) + qJ(3);
	t1 = [0, t50, t51, t50 * qJ(2) + t49 * t51 + 0; 0, -t51, t50, -t51 * qJ(2) + t49 * t50 + 0; 1, 0, 0, pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:31:41
	% EndTime: 2020-11-04 21:31:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (13->11), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->7)
	t57 = cos(qJ(1));
	t56 = cos(qJ(4));
	t55 = sin(qJ(1));
	t54 = sin(qJ(4));
	t53 = pkin(1) + qJ(3);
	t52 = pkin(7) - qJ(2);
	t1 = [t57 * t54, t57 * t56, -t55, -t52 * t55 + t53 * t57 + 0; t55 * t54, t55 * t56, t57, t52 * t57 + t53 * t55 + 0; t56, -t54, 0, pkin(3) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:31:41
	% EndTime: 2020-11-04 21:31:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (24->14), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t64 = cos(qJ(1));
	t63 = sin(qJ(1));
	t62 = qJ(4) + qJ(5);
	t61 = pkin(7) + pkin(8) - qJ(2);
	t60 = cos(t62);
	t59 = sin(t62);
	t58 = sin(qJ(4)) * pkin(4) + pkin(1) + qJ(3);
	t1 = [t64 * t59, t64 * t60, -t63, t58 * t64 - t61 * t63 + 0; t63 * t59, t63 * t60, t64, t58 * t63 + t61 * t64 + 0; t60, -t59, 0, cos(qJ(4)) * pkin(4) + pkin(3) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:31:41
	% EndTime: 2020-11-04 21:31:41
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (44->23), mult. (33->22), div. (0->0), fcn. (46->8), ass. (0->14)
	t70 = sin(qJ(6));
	t71 = sin(qJ(1));
	t78 = t71 * t70;
	t72 = cos(qJ(6));
	t77 = t71 * t72;
	t73 = cos(qJ(1));
	t76 = t73 * t70;
	t75 = t73 * t72;
	t69 = qJ(4) + qJ(5);
	t66 = sin(t69);
	t67 = cos(t69);
	t74 = pkin(5) * t66 - pkin(9) * t67 + sin(qJ(4)) * pkin(4) + pkin(1) + qJ(3);
	t68 = pkin(7) + pkin(8) - qJ(2);
	t1 = [t66 * t75 - t78, -t66 * t76 - t77, -t73 * t67, -t68 * t71 + t73 * t74 + 0; t66 * t77 + t76, -t66 * t78 + t75, -t71 * t67, t68 * t73 + t71 * t74 + 0; t67 * t72, -t67 * t70, t66, t67 * pkin(5) + t66 * pkin(9) + cos(qJ(4)) * pkin(4) + pkin(3) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end