% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRP3_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_jacobiR_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:44:46
% EndTime: 2019-02-26 20:44:46
% DurationCPUTime: 0.03s
% Computational Cost: add. (24->9), mult. (26->20), div. (0->0), fcn. (45->6), ass. (0->12)
t52 = sin(pkin(10));
t54 = sin(qJ(3));
t59 = t52 * t54;
t53 = cos(pkin(10));
t58 = t53 * t54;
t55 = cos(qJ(3));
t57 = t55 * t52;
t56 = t55 * t53;
t51 = qJ(1) + pkin(9);
t50 = cos(t51);
t49 = sin(t51);
t1 = [-t49 * t56 + t50 * t52, 0, -t50 * t58, 0, 0, 0; t49 * t52 + t50 * t56, 0, -t49 * t58, 0, 0, 0; 0, 0, t56, 0, 0, 0; t49 * t57 + t50 * t53, 0, t50 * t59, 0, 0, 0; t49 * t53 - t50 * t57, 0, t49 * t59, 0, 0, 0; 0, 0, -t57, 0, 0, 0; -t49 * t54, 0, t50 * t55, 0, 0, 0; t50 * t54, 0, t49 * t55, 0, 0, 0; 0, 0, t54, 0, 0, 0;];
JR_rot  = t1;
