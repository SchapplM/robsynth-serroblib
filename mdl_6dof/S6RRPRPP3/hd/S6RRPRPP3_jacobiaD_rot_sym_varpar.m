% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPP3
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6RRPRPP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPP3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP3_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPP3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
%% Function calls
if link_index == 0
	JaD_rot=S6RRPRPP3_jacobiaD_rot_0_sym_varpar(qJ, qJD, pkin);
elseif link_index == 1
	JaD_rot=S6RRPRPP3_jacobiaD_rot_1_sym_varpar(qJ, qJD, pkin);
elseif link_index == 2
	JaD_rot=S6RRPRPP3_jacobiaD_rot_2_sym_varpar(qJ, qJD, pkin);
elseif link_index == 3
	JaD_rot=S6RRPRPP3_jacobiaD_rot_3_sym_varpar(qJ, qJD, pkin);
elseif link_index == 4
	JaD_rot=S6RRPRPP3_jacobiaD_rot_4_sym_varpar(qJ, qJD, pkin);
elseif link_index == 5
	JaD_rot=S6RRPRPP3_jacobiaD_rot_5_sym_varpar(qJ, qJD, pkin);
elseif link_index == 6
	JaD_rot=S6RRPRPP3_jacobiaD_rot_6_sym_varpar(qJ, qJD, pkin);
else
	JaD_rot=NaN(3,6);
end