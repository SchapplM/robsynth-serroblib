% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPPRR5
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR5_jacobiR_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_jacobiR_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_jacobiR_rot_3_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:30:57
% EndTime: 2019-02-26 21:30:58
% DurationCPUTime: 0.02s
% Computational Cost: add. (8->6), mult. (28->12), div. (0->0), fcn. (48->6), ass. (0->15)
t61 = sin(qJ(2));
t62 = sin(qJ(1));
t68 = t62 * t61;
t63 = cos(qJ(2));
t67 = t62 * t63;
t64 = cos(qJ(1));
t66 = t64 * t61;
t65 = t64 * t63;
t60 = cos(pkin(6));
t59 = sin(pkin(6));
t58 = -t60 * t68 + t65;
t57 = t60 * t67 + t66;
t56 = t60 * t66 + t67;
t55 = t60 * t65 - t68;
t1 = [-t56, -t57, 0, 0, 0, 0; t58, t55, 0, 0, 0, 0; 0, t59 * t63, 0, 0, 0, 0; t64 * t59, 0, 0, 0, 0, 0; t62 * t59, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t55, t58, 0, 0, 0, 0; t57, t56, 0, 0, 0, 0; 0, t59 * t61, 0, 0, 0, 0;];
JR_rot  = t1;
