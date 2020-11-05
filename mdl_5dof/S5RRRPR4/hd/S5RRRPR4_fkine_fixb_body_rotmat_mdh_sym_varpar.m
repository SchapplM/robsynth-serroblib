% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRPR4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:42
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRRPR4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:42:37
	% EndTime: 2020-11-04 20:42:37
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:42:37
	% EndTime: 2020-11-04 20:42:37
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t52 = cos(qJ(1));
	t51 = sin(qJ(1));
	t1 = [t52, -t51, 0, 0; t51, t52, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:42:37
	% EndTime: 2020-11-04 20:42:37
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t55 = qJ(1) + qJ(2);
	t54 = cos(t55);
	t53 = sin(t55);
	t1 = [t54, -t53, 0, cos(qJ(1)) * pkin(1) + 0; t53, t54, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, pkin(6) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:42:37
	% EndTime: 2020-11-04 20:42:37
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t60 = cos(qJ(3));
	t59 = sin(qJ(3));
	t58 = qJ(1) + qJ(2);
	t57 = cos(t58);
	t56 = sin(t58);
	t1 = [t57 * t60, -t57 * t59, t56, t57 * pkin(2) + t56 * pkin(7) + cos(qJ(1)) * pkin(1) + 0; t56 * t60, -t56 * t59, -t57, t56 * pkin(2) - t57 * pkin(7) + sin(qJ(1)) * pkin(1) + 0; t59, t60, 0, pkin(6) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:42:37
	% EndTime: 2020-11-04 20:42:37
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (30->15), mult. (20->14), div. (0->0), fcn. (28->6), ass. (0->7)
	t64 = sin(qJ(3));
	t65 = cos(qJ(3));
	t66 = pkin(3) * t65 + qJ(4) * t64 + pkin(2);
	t63 = qJ(1) + qJ(2);
	t62 = cos(t63);
	t61 = sin(t63);
	t1 = [t62 * t65, t61, t62 * t64, cos(qJ(1)) * pkin(1) + t61 * pkin(7) + 0 + t66 * t62; t61 * t65, -t62, t61 * t64, sin(qJ(1)) * pkin(1) - t62 * pkin(7) + 0 + t66 * t61; t64, 0, -t65, t64 * pkin(3) - t65 * qJ(4) + pkin(5) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:42:37
	% EndTime: 2020-11-04 20:42:37
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (58->28), mult. (46->26), div. (0->0), fcn. (46->12), ass. (0->18)
	t79 = qJ(1) + qJ(2);
	t77 = qJ(3) + t79;
	t87 = sin(t77) / 0.2e1;
	t78 = -qJ(3) + t79;
	t86 = cos(t78) / 0.2e1;
	t85 = cos(qJ(5));
	t84 = sin(qJ(5));
	t83 = pkin(3) + pkin(4);
	t82 = pkin(7) - pkin(8);
	t81 = cos(qJ(3));
	t80 = sin(qJ(3));
	t76 = cos(t79);
	t75 = sin(t79);
	t73 = cos(t77);
	t72 = sin(t78);
	t68 = t80 * t85 - t81 * t84;
	t67 = -t80 * t84 - t81 * t85;
	t1 = [-t76 * t67, t76 * t68, -t75, t82 * t75 + cos(qJ(1)) * pkin(1) + t76 * pkin(2) + 0 + (t86 + t73 / 0.2e1) * t83 + (-t72 / 0.2e1 + t87) * qJ(4); -t75 * t67, t75 * t68, t76, -t82 * t76 + sin(qJ(1)) * pkin(1) + t75 * pkin(2) + 0 + (t72 / 0.2e1 + t87) * t83 + (t86 - t73 / 0.2e1) * qJ(4); t68, t67, 0, -t81 * qJ(4) + t83 * t80 + pkin(5) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end