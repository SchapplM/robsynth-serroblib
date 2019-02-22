% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR2
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:13
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR2_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_jacobiR_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:13:13
% EndTime: 2019-02-22 11:13:13
% DurationCPUTime: 0.06s
% Computational Cost: add. (23->9), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->12)
t57 = sin(pkin(11));
t59 = sin(qJ(1));
t64 = t59 * t57;
t58 = cos(pkin(11));
t63 = t59 * t58;
t60 = cos(qJ(1));
t62 = t60 * t57;
t61 = t60 * t58;
t56 = qJ(2) + pkin(10);
t55 = cos(t56);
t54 = sin(t56);
t1 = [-t55 * t63 + t62, -t54 * t61, 0, 0, 0, 0; t55 * t61 + t64, -t54 * t63, 0, 0, 0, 0; 0, t55 * t58, 0, 0, 0, 0; t55 * t64 + t61, t54 * t62, 0, 0, 0, 0; -t55 * t62 + t63, t54 * t64, 0, 0, 0, 0; 0, -t55 * t57, 0, 0, 0, 0; -t59 * t54, t60 * t55, 0, 0, 0, 0; t60 * t54, t59 * t55, 0, 0, 0, 0; 0, t54, 0, 0, 0, 0;];
JR_rot  = t1;
