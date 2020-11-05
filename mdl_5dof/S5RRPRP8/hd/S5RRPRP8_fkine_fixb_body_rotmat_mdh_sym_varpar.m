% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRP8 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:34
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRPRP8_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRP8_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:34:24
	% EndTime: 2020-11-04 20:34:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:34:24
	% EndTime: 2020-11-04 20:34:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t41 = cos(qJ(1));
	t40 = sin(qJ(1));
	t1 = [t41, -t40, 0, 0; t40, t41, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:34:24
	% EndTime: 2020-11-04 20:34:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t45 = cos(qJ(1));
	t44 = cos(qJ(2));
	t43 = sin(qJ(1));
	t42 = sin(qJ(2));
	t1 = [t45 * t44, -t45 * t42, t43, t45 * pkin(1) + t43 * pkin(6) + 0; t43 * t44, -t43 * t42, -t45, t43 * pkin(1) - t45 * pkin(6) + 0; t42, t44, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:34:24
	% EndTime: 2020-11-04 20:34:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (13->11), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->6)
	t50 = cos(qJ(1));
	t49 = cos(qJ(2));
	t48 = sin(qJ(1));
	t47 = sin(qJ(2));
	t46 = t49 * pkin(2) + t47 * qJ(3) + pkin(1);
	t1 = [t50 * t49, t48, t50 * t47, t48 * pkin(6) + t46 * t50 + 0; t48 * t49, -t50, t48 * t47, -t50 * pkin(6) + t46 * t48 + 0; t47, 0, -t49, t47 * pkin(2) - t49 * qJ(3) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:34:24
	% EndTime: 2020-11-04 20:34:24
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (23->16), mult. (26->16), div. (0->0), fcn. (40->6), ass. (0->12)
	t63 = cos(qJ(4));
	t62 = sin(qJ(4));
	t61 = pkin(2) + pkin(3);
	t60 = pkin(6) - pkin(7);
	t59 = cos(qJ(1));
	t58 = cos(qJ(2));
	t57 = sin(qJ(1));
	t56 = sin(qJ(2));
	t53 = t56 * qJ(3) + t61 * t58 + pkin(1);
	t52 = t56 * t63 - t58 * t62;
	t51 = -t56 * t62 - t58 * t63;
	t1 = [-t59 * t51, t59 * t52, -t57, t53 * t59 + t60 * t57 + 0; -t57 * t51, t57 * t52, t59, t53 * t57 - t60 * t59 + 0; t52, t51, 0, -t58 * qJ(3) + t61 * t56 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:34:24
	% EndTime: 2020-11-04 20:34:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (31->18), mult. (32->18), div. (0->0), fcn. (46->6), ass. (0->13)
	t77 = cos(qJ(1));
	t76 = cos(qJ(2));
	t75 = cos(qJ(4));
	t74 = sin(qJ(1));
	t73 = sin(qJ(2));
	t72 = sin(qJ(4));
	t71 = qJ(5) - pkin(6) + pkin(7);
	t68 = t72 * pkin(4) + qJ(3);
	t67 = t75 * pkin(4) + pkin(2) + pkin(3);
	t66 = -t76 * t72 + t73 * t75;
	t65 = -t73 * t72 - t76 * t75;
	t64 = t67 * t76 + t68 * t73 + pkin(1);
	t1 = [-t77 * t65, t77 * t66, -t74, t64 * t77 - t71 * t74 + 0; -t74 * t65, t74 * t66, t77, t64 * t74 + t71 * t77 + 0; t66, t65, 0, t67 * t73 - t68 * t76 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end