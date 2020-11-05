% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRPPR4 (for one body)
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
% Datum: 2020-11-04 19:57
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PRPPR4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:57:14
	% EndTime: 2020-11-04 19:57:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:57:14
	% EndTime: 2020-11-04 19:57:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t52 = cos(pkin(7));
	t51 = sin(pkin(7));
	t1 = [t52, -t51, 0, 0; t51, t52, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:57:14
	% EndTime: 2020-11-04 19:57:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t55 = pkin(7) + qJ(2);
	t54 = cos(t55);
	t53 = sin(t55);
	t1 = [t54, -t53, 0, cos(pkin(7)) * pkin(1) + 0; t53, t54, 0, sin(pkin(7)) * pkin(1) + 0; 0, 0, 1, pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:57:14
	% EndTime: 2020-11-04 19:57:15
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t60 = cos(pkin(8));
	t59 = sin(pkin(8));
	t58 = pkin(7) + qJ(2);
	t57 = cos(t58);
	t56 = sin(t58);
	t1 = [t57 * t60, -t57 * t59, t56, t57 * pkin(2) + t56 * qJ(3) + cos(pkin(7)) * pkin(1) + 0; t56 * t60, -t56 * t59, -t57, t56 * pkin(2) - t57 * qJ(3) + sin(pkin(7)) * pkin(1) + 0; t59, t60, 0, pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:57:15
	% EndTime: 2020-11-04 19:57:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (30->15), mult. (20->14), div. (0->0), fcn. (28->6), ass. (0->7)
	t64 = sin(pkin(8));
	t65 = cos(pkin(8));
	t66 = pkin(3) * t65 + qJ(4) * t64 + pkin(2);
	t63 = pkin(7) + qJ(2);
	t62 = cos(t63);
	t61 = sin(t63);
	t1 = [t62 * t65, t61, t62 * t64, cos(pkin(7)) * pkin(1) + t61 * qJ(3) + 0 + t66 * t62; t61 * t65, -t62, t61 * t64, sin(pkin(7)) * pkin(1) - t62 * qJ(3) + 0 + t66 * t61; t64, 0, -t65, t64 * pkin(3) - t65 * qJ(4) + pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:57:15
	% EndTime: 2020-11-04 19:57:15
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (58->28), mult. (46->26), div. (0->0), fcn. (46->12), ass. (0->18)
	t79 = pkin(7) + qJ(2);
	t75 = pkin(8) + t79;
	t87 = sin(t75) / 0.2e1;
	t76 = -pkin(8) + t79;
	t86 = cos(t76) / 0.2e1;
	t85 = cos(qJ(5));
	t84 = sin(qJ(5));
	t83 = pkin(3) + pkin(4);
	t82 = qJ(3) - pkin(6);
	t81 = cos(pkin(8));
	t80 = sin(pkin(8));
	t78 = cos(t79);
	t77 = sin(t79);
	t71 = cos(t75);
	t70 = sin(t76);
	t68 = t80 * t85 - t81 * t84;
	t67 = -t80 * t84 - t81 * t85;
	t1 = [-t78 * t67, t78 * t68, -t77, t82 * t77 + cos(pkin(7)) * pkin(1) + t78 * pkin(2) + 0 + (t86 + t71 / 0.2e1) * t83 + (-t70 / 0.2e1 + t87) * qJ(4); -t77 * t67, t77 * t68, t78, -t82 * t78 + sin(pkin(7)) * pkin(1) + t77 * pkin(2) + 0 + (t70 / 0.2e1 + t87) * t83 + (t86 - t71 / 0.2e1) * qJ(4); t68, t67, 0, -t81 * qJ(4) + t83 * t80 + pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end